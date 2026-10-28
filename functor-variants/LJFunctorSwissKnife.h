/**
 * @file LJFunctorSwissKnife.h
 * @brief High-performance Lennard-Jones Functor utilizing multi-regime mixing:
 *        1. Homogeneous: Pure scalar broadcast for single-component or pure cells.
 *        2. Row Hoist: Direct table row access for small cells (N <= 32).
 *        3. Compact Preload Cache: Aligned SIMD vector loads for large cells (N > 32).
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
class LJFunctorSwissKnife
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorSwissKnife<Particle_T, applyShift, useMixing, useNewton3,
                                                                      calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  enum LJMixParam : size_t { SigmaSquared = 0, Epsilon24 = 1, Shift6 = 2 };
  using PreloadCacheType = SoAMixingPreloadCache<SoAFloatPrecision, applyShift ? 3 : 2, 4, 256>;

  static constexpr size_t SmallCellThreshold = 32;

 public:
  LJFunctorSwissKnife() = delete;

 private:
  explicit LJFunctorSwissKnife(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorSwissKnife<Particle_T, applyShift, useMixing, useNewton3,
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
  explicit LJFunctorSwissKnife(double cutoff) : LJFunctorSwissKnife(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  explicit LJFunctorSwissKnife(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorSwissKnife(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    initMixingTableFromPPL(particlePropertiesLibrary);
  }

  void initMixingTableFromPPL(const ParticlePropertiesLibrary<double, size_t> &ppl) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _mixingTable.build(nTypes, [&](size_t i, size_t j) -> std::array<SoAFloatPrecision, 3> {
      const auto data = ppl.getLJMixingData(i, j);
      return {static_cast<SoAFloatPrecision>(data.sigmaSquared),
              static_cast<SoAFloatPrecision>(data.epsilon24),
              static_cast<SoAFloatPrecision>(data.shift6)};
    });
    if (nTypes > 0) {
      _sigmaSquared = _mixingTable.template get<SigmaSquared>(0, 0);
      _epsilon24 = _mixingTable.template get<Epsilon24>(0, 0);
      if constexpr (applyShift) {
        _shift6 = _mixingTable.template get<Shift6>(0, 0);
      }
    }
  }

  std::string getName() final { return "LJFunctorSwissKnife"; }
  bool isRelevantForTuning() final { return relevantForTuning; }
  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }
  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

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

    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24;
    auto shift6 = _shift6;
    if constexpr (useMixing) {
      const size_t tI = i.getTypeId();
      const size_t tJ = j.getTypeId();
      sigmaSquared = _mixingTable.template get<LJMixParam::SigmaSquared>(tI, tJ);
      epsilon24 = _mixingTable.template get<LJMixParam::Epsilon24>(tI, tJ);
      if constexpr (applyShift) {
        shift6 = _mixingTable.template get<LJMixParam::Shift6>(tI, tJ);
      }
    }
    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
    if (dr2 > _cutoffSquared) return;

    double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac;
    i.addF(f);
    if (newton3) j.subF(f);

    if constexpr (calculateGlobals) {
      const auto threadnum = autopas::autopas_get_thread_num();
      auto virial = dr * f;
      double potentialEnergy6 = epsilon24 * lj12m6 + shift6;
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
      }
    }
  }

  // --- Parameter Accessors for unified kernel templating ---

  struct HomogeneousAccessor {
    SoAFloatPrecision sigmaSquared, epsilon24, shift6;
    [[gnu::always_inline]] void fetch(size_t /*j*/, SoAFloatPrecision &s, SoAFloatPrecision &e,
                                      SoAFloatPrecision &sh) const noexcept {
      s = sigmaSquared;
      e = epsilon24;
      sh = shift6;
    }
  };

  struct RowHoistAccessor {
    const SoAFloatPrecision *const __restrict__ sigmaRow;
    const SoAFloatPrecision *const __restrict__ epsRow;
    const SoAFloatPrecision *const __restrict__ shiftRow;
    const size_t *const __restrict__ typeptrJ;
    [[gnu::always_inline]] void fetch(size_t j, SoAFloatPrecision &s, SoAFloatPrecision &e,
                                      SoAFloatPrecision &sh) const noexcept {
      const size_t tJ = typeptrJ[j];
      s = sigmaRow[tJ];
      e = epsRow[tJ];
      if constexpr (applyShift) sh = shiftRow[tJ];
    }
  };

  struct PreloadAccessor {
    const SoAFloatPrecision *const __restrict__ sigmaArr;
    const SoAFloatPrecision *const __restrict__ epsArr;
    const SoAFloatPrecision *const __restrict__ shiftArr;
    [[gnu::always_inline]] void fetch(size_t j, SoAFloatPrecision &s, SoAFloatPrecision &e,
                                      SoAFloatPrecision &sh) const noexcept {
      s = sigmaArr[j];
      e = epsArr[j];
      if constexpr (applyShift) sh = shiftArr[j];
    }
  };

  // --- SoAFunctorSingle ---

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    if (soa.size() == 0) return;

    if constexpr (useMixing) {
      if (_mixingTable.numTypes <= 1 or isCellHomogeneous(soa)) {
        const auto t0 = soa.template begin<Particle_T::AttributeNames::typeId>()[0];
        HomogeneousAccessor acc{
            _mixingTable.template get<SigmaSquared>(t0, t0),
            _mixingTable.template get<Epsilon24>(t0, t0),
            applyShift ? _mixingTable.template get<Shift6>(t0, t0) : SoAFloatPrecision{0.}};
        SoAFunctorSingleLoop(soa, newton3, [&](size_t /*i*/) { return acc; });
      } else if (soa.size() <= SmallCellThreshold) {
        // Small cell regime: Direct Row Hoisting (0 cache overhead)
        const auto *const __restrict__ typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
        SoAFunctorSingleLoop(soa, newton3, [&](size_t i) {
          const size_t tI = typeptr[i];
          return RowHoistAccessor{
              _mixingTable.template row<SigmaSquared>(tI),
              _mixingTable.template row<Epsilon24>(tI),
              applyShift ? _mixingTable.template row<Shift6>(tI) : nullptr,
              typeptr};
        });
      } else {
        // Large cell regime: Compact Preload Cache
        const auto *const __restrict__ typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
        PreloadCacheType cache{_mixingTable.numTypes, soa.size()};
        SoAFunctorSingleLoop(soa, newton3, [&](size_t i) {
          const size_t tI = typeptr[i];
          if constexpr (applyShift) {
            auto [s, e, sh] = cache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                _mixingTable, tI, typeptr, soa.size());
            return PreloadAccessor{s, e, sh};
          } else {
            auto [s, e] = cache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                _mixingTable, tI, typeptr, soa.size());
            return PreloadAccessor{s, e, nullptr};
          }
        });
      }
    } else {
      HomogeneousAccessor acc{_sigmaSquared, _epsilon24, _shift6};
      SoAFunctorSingleLoop(soa, newton3, [&](size_t /*i*/) { return acc; });
    }
  }

  template <typename AccessorFactory>
  void SoAFunctorSingleLoop(autopas::SoAView<SoAArraysType> soa, bool newton3, AccessorFactory &&getAccessor) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    for (unsigned int i = 0; i < soa.size(); ++i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) continue;

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      const auto accessor = getAccessor(i);

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsSum)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        SoAFloatPrecision sigmaSquared, epsilon24, shift6;
        accessor.fetch(j, sigmaSquared, epsilon24, shift6);

        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;
        const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
        const SoAFloatPrecision lj4 = lj2 * lj2;
        const SoAFloatPrecision lj6 = lj4 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision lj12m6alj12 = lj12m6 + lj12;
        const SoAFloatPrecision lj12m6alj12e = lj12m6alj12 * epsilon24;
        const SoAFloatPrecision fac = mask ? lj12m6alj12e * invdr2 : 0.;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if constexpr (calculateGlobals) {
          if (mask) {
            const SoAFloatPrecision potentialEnergy6 = epsilon24 * lj12m6 + shift6;
            const bool ownedI = ownedStatePtr[i] == autopas::OwnershipState::owned;
            const bool ownedJ = ownedStateJ == autopas::OwnershipState::owned;
            if (ownedI) {
              potentialEnergySum += potentialEnergy6;
              virialSumX += drx * fx;
              virialSumY += dry * fy;
              virialSumZ += drz * fz;
            }
            if (newton3 and ownedJ) {
              potentialEnergySum += potentialEnergy6;
              virialSumX += drx * fx;
              virialSumY += dry * fy;
              virialSumZ += drz * fz;
            }
          }
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  // --- SoAFunctorPair ---

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) return;

    if constexpr (useMixing) {
      const auto *const __restrict__ typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
      const auto *const __restrict__ typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

      if (_mixingTable.numTypes <= 1 or (isCellHomogeneous(soa1) and isCellHomogeneous(soa2) and typeptr1[0] == typeptr2[0])) {
        const auto t0 = typeptr1[0];
        HomogeneousAccessor acc{
            _mixingTable.template get<SigmaSquared>(t0, t0),
            _mixingTable.template get<Epsilon24>(t0, t0),
            applyShift ? _mixingTable.template get<Shift6>(t0, t0) : SoAFloatPrecision{0.}};
        SoAFunctorPairLoop(soa1, soa2, newton3, [&](size_t /*i*/) { return acc; });
      } else if (soa2.size() <= SmallCellThreshold) {
        // Small cell regime: Row Hoisting
        SoAFunctorPairLoop(soa1, soa2, newton3, [&](size_t i) {
          const size_t tI = typeptr1[i];
          return RowHoistAccessor{
              _mixingTable.template row<SigmaSquared>(tI),
              _mixingTable.template row<Epsilon24>(tI),
              applyShift ? _mixingTable.template row<Shift6>(tI) : nullptr,
              typeptr2};
        });
      } else {
        // Large cell regime: Compact Preload Cache
        PreloadCacheType cache{_mixingTable.numTypes, soa2.size()};
        SoAFunctorPairLoop(soa1, soa2, newton3, [&](size_t i) {
          const size_t tI = typeptr1[i];
          if constexpr (applyShift) {
            auto [s, e, sh] = cache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                _mixingTable, tI, typeptr2, soa2.size());
            return PreloadAccessor{s, e, sh};
          } else {
            auto [s, e] = cache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                _mixingTable, tI, typeptr2, soa2.size());
            return PreloadAccessor{s, e, nullptr};
          }
        });
      }
    } else {
      HomogeneousAccessor acc{_sigmaSquared, _epsilon24, _shift6};
      SoAFunctorPairLoop(soa1, soa2, newton3, [&](size_t /*i*/) { return acc; });
    }
  }

  template <typename AccessorFactory>
  void SoAFunctorPairLoop(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                          bool newton3, AccessorFactory &&getAccessor) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fx1ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy1ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz1ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    SoAFloatPrecision *const __restrict fx2ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy2ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz2ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) continue;

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      const auto accessor = getAccessor(i);

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsSum)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        SoAFloatPrecision sigmaSquared, epsilon24, shift6;
        accessor.fetch(j, sigmaSquared, epsilon24, shift6);

        const auto ownedStateJ = ownedStatePtr2[j];

        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;
        const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
        const SoAFloatPrecision lj4 = lj2 * lj2;
        const SoAFloatPrecision lj6 = lj4 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision lj12m6alj12 = lj12m6 + lj12;
        const SoAFloatPrecision lj12m6alj12e = lj12m6alj12 * epsilon24;
        const SoAFloatPrecision fac = mask ? lj12m6alj12e * invdr2 : 0.;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }

        if constexpr (calculateGlobals) {
          if (mask) {
            const SoAFloatPrecision potentialEnergy6 = epsilon24 * lj12m6 + shift6;
            const bool ownedI = ownedStatePtr1[i] == autopas::OwnershipState::owned;
            const bool ownedJ = ownedStateJ == autopas::OwnershipState::owned;
            if (ownedI) {
              potentialEnergySum += potentialEnergy6;
              virialSumX += drx * fx;
              virialSumY += dry * fy;
              virialSumZ += drz * fz;
            }
            if (newton3 and ownedJ) {
              potentialEnergySum += potentialEnergy6;
              virialSumX += drx * fx;
              virialSumY += dry * fy;
              virialSumZ += drz * fz;
            }
          }
        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        std::span<const size_t> neighborList, bool newton3) final {
    // Verlet fallback implementation using direct table access
    if (soa.size() == 0 or neighborList.empty()) return;
    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) return;

    SoAFloatPrecision fxacc = 0, fyacc = 0, fzacc = 0;
    const size_t typeI = useMixing ? typeptr[indexFirst] : 0;
    const auto *const __restrict sigmaRow = useMixing ? _mixingTable.template row<SigmaSquared>(typeI) : nullptr;
    const auto *const __restrict epsRow = useMixing ? _mixingTable.template row<Epsilon24>(typeI) : nullptr;
    const auto *const __restrict shiftRow = (useMixing and applyShift) ? _mixingTable.template row<Shift6>(typeI) : nullptr;

    const size_t numNeighbors = neighborList.size();
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
    for (size_t n = 0; n < numNeighbors; ++n) {
      const size_t j = neighborList[n];
      const auto ownedStateJ = ownedStatePtr[j];
      SoAFloatPrecision sigmaSquared = _sigmaSquared;
      SoAFloatPrecision epsilon24 = _epsilon24;
      SoAFloatPrecision shift6 = _shift6;
      if constexpr (useMixing) {
        const size_t typeJ = typeptr[j];
        sigmaSquared = sigmaRow[typeJ];
        epsilon24 = epsRow[typeJ];
        if constexpr (applyShift) shift6 = shiftRow[typeJ];
      }

      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];
      const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;
      const bool mask = dr2 <= _cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;
      const SoAFloatPrecision invdr2 = 1. / dr2;
      const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
      const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
      const SoAFloatPrecision lj12 = lj6 * lj6;
      const SoAFloatPrecision fac = mask ? epsilon24 * (lj12 + lj12 - lj6) * invdr2 * invdr2 : 0.;

      const SoAFloatPrecision fx = drx * fac;
      const SoAFloatPrecision fy = dry * fac;
      const SoAFloatPrecision fz = drz * fac;

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;

      if (newton3) {
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }
    }

    fxptr[indexFirst] += fxacc;
    fyptr[indexFirst] += fyacc;
    fzptr[indexFirst] += fzacc;
  }

  // --- Globals & Traversal Boilerplate ---

  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __pad{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }
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
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    if constexpr (calculateGlobals) {
      for (auto &data : _aosThreadDataGlobals) data.setZero();
    }
  }

  void endTraversal(bool /*newton3*/) final {
    if (_postProcessed) return;
    if constexpr (calculateGlobals) {
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum[0] += data.virialSum[0];
        _virialSum[1] += data.virialSum[1];
        _virialSum[2] += data.virialSum[2];
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

  [[gnu::always_inline]] static bool isCellHomogeneous(const autopas::SoAView<SoAArraysType> &soa) noexcept {
    const auto numParticles = soa.size();
    if (numParticles <= 1) return true;
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
    const auto firstType = typeptr[0];
    for (size_t i = 1; i < numParticles; ++i) {
      if (typeptr[i] != firstType) return false;
    }
    return true;
  }

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
