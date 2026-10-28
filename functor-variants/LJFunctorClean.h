/**
 * @file LJFunctorClean.h
 * @brief Clean, maintainable Lennard-Jones Functor demonstrating how modern AutoPas functors
 *        can remain concise, readable, and under 350 lines with zero mixing boilerplate.
 */

#pragma once

#include <cstddef>
#include <span>
#include <tuple>
#include "MixingTable.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorClean
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorClean<Particle_T, applyShift, useMixing, useNewton3,
                                                                 calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  enum LJMixParam : size_t { SigmaSquared = 0, Epsilon24 = 1, Shift6 = 2 };
  using CacheType = SoAMixingPreloadCache<SoAFloatPrecision, applyShift ? 3 : 2, 4, 256>;

 public:
  LJFunctorClean() = delete;

 private:
  explicit LJFunctorClean(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorClean<Particle_T, applyShift, useMixing, useNewton3,
                                                            calculateGlobals, countFLOPs, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  explicit LJFunctorClean(double cutoff) : LJFunctorClean(cutoff, nullptr) {
    static_assert(not useMixing, "Mixing without ParticlePropertiesLibrary is not permitted.");
  }

  explicit LJFunctorClean(double cutoff, ParticlePropertiesLibrary<double, size_t> &ppl)
      : LJFunctorClean(cutoff, nullptr) {
    static_assert(useMixing, "Cannot use ParticlePropertiesLibrary with mixing disabled.");
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _mixingTable.build(nTypes, [&](size_t i, size_t j) -> std::array<SoAFloatPrecision, 3> {
      const auto d = ppl.getLJMixingData(i, j);
      return {static_cast<SoAFloatPrecision>(d.sigmaSquared),
              static_cast<SoAFloatPrecision>(d.epsilon24),
              static_cast<SoAFloatPrecision>(d.shift6)};
    });
    if (nTypes > 0) {
      _sigmaSquared = _mixingTable.template get<SigmaSquared>(0, 0);
      _epsilon24 = _mixingTable.template get<Epsilon24>(0, 0);
      if constexpr (applyShift) _shift6 = _mixingTable.template get<Shift6>(0, 0);
    }
  }

  std::string getName() final { return "LJFunctorClean"; }
  bool isRelevantForTuning() final { return relevantForTuning; }
  bool allowsNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Off; }
  bool allowsNonNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Only; }

  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getNeededAttr(std::false_type) { return getNeededAttr(); }

  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) return;

    auto s2 = _sigmaSquared, e24 = _epsilon24, sh6 = _shift6;
    if constexpr (useMixing) {
      const size_t tI = i.getTypeId(), tJ = j.getTypeId();
      s2 = _mixingTable.template get<SigmaSquared>(tI, tJ);
      e24 = _mixingTable.template get<Epsilon24>(tI, tJ);
      if constexpr (applyShift) sh6 = _mixingTable.template get<Shift6>(tI, tJ);
    }
    const auto dr = i.getR() - j.getR();
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
    if (dr2 > _cutoffSquared) return;

    const double invdr2 = 1. / dr2;
    const double lj6 = (s2 * invdr2) * (s2 * invdr2) * (s2 * invdr2);
    const double lj12 = lj6 * lj6;
    const double fac = e24 * (lj12 + lj12 - lj6) * invdr2;
    const auto f = dr * fac;
    i.addF(f);
    if (newton3) j.subF(f);

    if constexpr (calculateGlobals) {
      const auto threadnum = autopas::autopas_get_thread_num();
      const double pe = e24 * (lj12 - lj6) + sh6;
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += pe;
        _aosThreadDataGlobals[threadnum].virialSum += dr * f;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += pe;
        _aosThreadDataGlobals[threadnum].virialSum += dr * f;
      }
    }
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    if (soa.size() == 0) return;
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict x = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict owned = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz = soa.template begin<Particle_T::AttributeNames::forceZ>();

    CacheType cache{_mixingTable.numTypes, soa.size()};

    SoAFloatPrecision peSum = 0., virX = 0., virY = 0., virZ = 0.;
    const auto cutoffSq = _cutoffSquared;

    for (unsigned int i = 0; i < soa.size(); ++i) {
      if (owned[i] == autopas::OwnershipState::dummy) continue;

      SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;
      const SoAFloatPrecision *sArr = nullptr, *eArr = nullptr, *shArr = nullptr;

      if constexpr (useMixing) {
        if constexpr (applyShift) {
          std::tie(sArr, eArr, shArr) = cache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
              _mixingTable, typeptr[i], typeptr, soa.size());
        } else {
          std::tie(sArr, eArr) = cache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
              _mixingTable, typeptr[i], typeptr, soa.size());
        }
      }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        const auto s2 = useMixing ? sArr[j] : _sigmaSquared;
        const auto e24 = useMixing ? eArr[j] : _epsilon24;
        const auto sh6 = (useMixing and applyShift) ? shArr[j] : _shift6;

        const auto drx = x[i] - x[j], dry = y[i] - y[j], drz = z[i] - z[j];
        const auto dr2 = drx * drx + dry * dry + drz * drz;
        const bool mask = dr2 <= cutoffSq and owned[j] != autopas::OwnershipState::dummy;

        const auto invdr2 = 1. / dr2;
        const auto lj2 = s2 * invdr2;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : 0.;

        const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
        fxacc += dfx; fyacc += dfy; fzacc += dfz;

        fx[j] -= dfx;
        fy[j] -= dfy;
        fz[j] -= dfz;
        if constexpr (calculateGlobals) {
          if (mask) {
            const auto pe = e24 * (lj12 - lj6) + sh6;
            if (owned[i] == autopas::OwnershipState::owned) {
              peSum += pe; virX += drx * dfx; virY += dry * dfy; virZ += drz * dfz;
            }
            if (newton3 and owned[j] == autopas::OwnershipState::owned) {
              peSum += pe; virX += drx * dfx; virY += dry * dfy; virZ += drz * dfz;
            }
          }
        }
      }

      fx[i] += fxacc; fy[i] += fyacc; fz[i] += fzacc;
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += peSum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virZ;
    }
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) return;
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict x1 = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1 = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1 = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict x2 = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2 = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2 = soa2.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict owned1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict owned2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict type1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    const auto *const __restrict type2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx1 = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1 = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1 = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fx2 = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2 = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2 = soa2.template begin<Particle_T::AttributeNames::forceZ>();

    CacheType cache{_mixingTable.numTypes, soa2.size()};

    SoAFloatPrecision peSum = 0., virX = 0., virY = 0., virZ = 0.;
    const auto cutoffSq = _cutoffSquared;

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (owned1[i] == autopas::OwnershipState::dummy) continue;

      SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;
      const SoAFloatPrecision *sArr = nullptr, *eArr = nullptr, *shArr = nullptr;

      if constexpr (useMixing) {
        if constexpr (applyShift) {
          std::tie(sArr, eArr, shArr) = cache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
              _mixingTable, type1[i], type2, soa2.size());
        } else {
          std::tie(sArr, eArr) = cache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
              _mixingTable, type1[i], type2, soa2.size());
        }
      }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto s2 = useMixing ? sArr[j] : _sigmaSquared;
        const auto e24 = useMixing ? eArr[j] : _epsilon24;
        const auto sh6 = (useMixing and applyShift) ? shArr[j] : _shift6;

        const auto drx = x1[i] - x2[j], dry = y1[i] - y2[j], drz = z1[i] - z2[j];
        const auto dr2 = drx * drx + dry * dry + drz * drz;
        const bool mask = dr2 <= cutoffSq and owned2[j] != autopas::OwnershipState::dummy;

        const auto invdr2 = 1. / dr2;
        const auto lj2 = s2 * invdr2;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : 0.;

        const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
        fxacc += dfx; fyacc += dfy; fzacc += dfz;

        if (newton3) {
          fx2[j] -= dfx; fy2[j] -= dfy; fz2[j] -= dfz;
        }
        if constexpr (calculateGlobals) {
          if (mask) {
            const auto pe = e24 * (lj12 - lj6) + sh6;
            if (owned1[i] == autopas::OwnershipState::owned) {
              peSum += pe; virX += drx * dfx; virY += dry * dfy; virZ += drz * dfz;
            }
            if (newton3 and owned2[j] == autopas::OwnershipState::owned) {
              peSum += pe; virX += drx * dfx; virY += dry * dfy; virZ += drz * dfz;
            }
          }
        }
      }

      fx1[i] += fxacc; fy1[i] += fyacc; fz1[i] += fzacc;
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += peSum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virZ;
    }
  }

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        std::span<const size_t> neighborList, bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    const auto *const __restrict x = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict owned = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
    auto *const __restrict fx = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz = soa.template begin<Particle_T::AttributeNames::forceZ>();

    if (owned[indexFirst] == autopas::OwnershipState::dummy) return;
    SoAFloatPrecision fxacc = 0, fyacc = 0, fzacc = 0;
    const size_t tI = useMixing ? typeptr[indexFirst] : 0;
    const auto *const __restrict sRow = useMixing ? _mixingTable.template row<SigmaSquared>(tI) : nullptr;
    const auto *const __restrict eRow = useMixing ? _mixingTable.template row<Epsilon24>(tI) : nullptr;
    const auto *const __restrict shRow = (useMixing and applyShift) ? _mixingTable.template row<Shift6>(tI) : nullptr;

    const size_t nNeighbors = neighborList.size();
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
    for (size_t n = 0; n < nNeighbors; ++n) {
      const size_t j = neighborList[n];
      const auto s2 = useMixing ? sRow[typeptr[j]] : _sigmaSquared;
      const auto e24 = useMixing ? eRow[typeptr[j]] : _epsilon24;
      const auto drx = x[indexFirst] - x[j], dry = y[indexFirst] - y[j], drz = z[indexFirst] - z[j];
      const auto dr2 = drx * drx + dry * dry + drz * drz;
      const bool mask = dr2 <= _cutoffSquared and owned[j] != autopas::OwnershipState::dummy;
      const auto invdr2 = 1. / dr2;
      const auto lj2 = s2 * invdr2;
      const auto lj6 = lj2 * lj2 * lj2;
      const auto lj12 = lj6 * lj6;
      const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 * invdr2 : 0.;
      const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
      fxacc += dfx; fyacc += dfy; fzacc += dfz;
      if (newton3) { fx[j] -= dfx; fy[j] -= dfy; fz[j] -= dfz; }
    }
    fx[indexFirst] += fxacc; fy[indexFirst] += fyacc; fz[indexFirst] += fzacc;
  }

  // --- Globals & Traversal Boilerplate ---
  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __pad{} {}
    void setZero() { virialSum = {0., 0., 0.}; potentialEnergySum = 0.; }
    std::array<double, 3> virialSum;
    double potentialEnergySum;
   private:
    double __pad[(64 - 4 * sizeof(double)) / sizeof(double)];
  };

  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __pad{} {}
    void setZero() {}
   private:
    double __pad[64 / sizeof(double)];
  };

  void initTraversal() final {
    _potentialEnergySum = 0.; _virialSum = {0., 0., 0.}; _postProcessed = false;
    if constexpr (calculateGlobals) for (auto &d : _aosThreadDataGlobals) d.setZero();
  }

  void endTraversal(bool /*newton3*/) final {
    if (_postProcessed) return;
    if constexpr (calculateGlobals) {
      for (const auto &d : _aosThreadDataGlobals) {
        _potentialEnergySum += d.potentialEnergySum;
        _virialSum[0] += d.virialSum[0];
        _virialSum[1] += d.virialSum[1];
        _virialSum[2] += d.virialSum[2];
      }
      _potentialEnergySum *= 0.5;
      _virialSum[0] *= 0.5;
      _virialSum[1] *= 0.5;
      _virialSum[2] *= 0.5;
    }
    _postProcessed = true;
  }

  double getPotentialEnergy() const { return _potentialEnergySum; }
  double getVirial() const { return _virialSum[0] + _virialSum[1] + _virialSum[2]; }

 private:
  CompactMixingTable<2, 3, SoAFloatPrecision, 8> _mixingTable;
  const double _cutoffSquared;
  SoAFloatPrecision _sigmaSquared{1.};
  SoAFloatPrecision _epsilon24{24.};
  SoAFloatPrecision _shift6{0.};
  double _potentialEnergySum;
  std::array<double, 3> _virialSum;
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};
  bool _postProcessed;
};

}  // namespace mdLib
