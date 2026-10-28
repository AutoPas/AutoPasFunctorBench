/**
 * @file LJFunctorNewMixing.h
 * @brief Simplified, zero-duplication LJFunctor using CompactMixingTable + SoAMixingPreloadCache.
 */

#pragma once

#include <tuple>
#include <span>

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
class LJFunctorNewMixing
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorNewMixing<Particle_T, applyShift, useMixing, useNewton3,
                                                                     calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  enum LJMixParam : size_t { SigmaSquared = 0, Epsilon24 = 1, Shift6 = 2 };
  using PreloadCacheType = SoAMixingPreloadCache<SoAFloatPrecision, applyShift ? 3 : 2>;

 public:
  LJFunctorNewMixing() = delete;

 private:
  explicit LJFunctorNewMixing(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorNewMixing<Particle_T, applyShift, useMixing, useNewton3,
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
  explicit LJFunctorNewMixing(double cutoff) : LJFunctorNewMixing(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  explicit LJFunctorNewMixing(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorNewMixing(cutoff, nullptr) {
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
    if (nTypes >= 1) {
      _sigmaSquared = _mixingTable.template get<LJMixParam::SigmaSquared>(0, 0);
      _epsilon24 = _mixingTable.template get<LJMixParam::Epsilon24>(0, 0);
      if constexpr (applyShift) {
        _shift6 = _mixingTable.template get<LJMixParam::Shift6>(0, 0);
      }
    }
  }

  std::string getName() final { return "LJFunctorNewMixing"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

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

    if (dr2 > _cutoffSquared) {
      return;
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac;
    i.addF(f);
    if (newton3) {
      j.subF(f);
    }

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGlobals) {
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
      if constexpr (countFLOPs) {
        if (newton3) {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
        } else {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
        }
      }
    }
  }

  template <typename SoAViewType>
  [[gnu::always_inline]] static bool isCellHomogeneous(const SoAViewType &soa) noexcept {
    const size_t sz = soa.size();
    if (sz <= 1) return true;
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
    const auto firstType = typeptr[0];
    for (size_t k = 1; k < sz; ++k) {
      if (typeptr[k] != firstType) return false;
    }
    return true;
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() == 0) return;
    if constexpr (useMixing) {
      const auto *const typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
      if (_mixingTable.numTypes <= 1 or isCellHomogeneous(soa)) {
        const auto t0 = typeptr[0];
        SoAFunctorSingleImpl<true>(soa, newton3,
                                   _mixingTable.template get<SigmaSquared>(t0, t0),
                                   _mixingTable.template get<Epsilon24>(t0, t0),
                                   applyShift ? _mixingTable.template get<Shift6>(t0, t0) : 0.);
      } else {
        SoAFunctorSingleImpl<false>(soa, newton3, 0., 0., 0.);
      }
    } else {
      SoAFunctorSingleImpl<true>(soa, newton3, _sigmaSquared, _epsilon24, _shift6);
    }
  }

  template <bool isHomogeneous>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa, bool newton3,
                            SoAFloatPrecision h_sigmaSquared,
                            SoAFloatPrecision h_epsilon24,
                            SoAFloatPrecision h_shift6) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    [[maybe_unused]] PreloadCacheType preloadCache;
    if constexpr (!isHomogeneous and useMixing) {
      preloadCache = PreloadCacheType{_mixingTable.numTypes, soa.size()};
    }

    for (unsigned int i = 0; i < soa.size(); ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      const SoAFloatPrecision *sigmaSquareds = nullptr;
      const SoAFloatPrecision *epsilon24s = nullptr;
      const SoAFloatPrecision *shift6s = nullptr;
      if constexpr (!isHomogeneous and useMixing) {
        if constexpr (applyShift) {
          std::tie(sigmaSquareds, epsilon24s, shift6s) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeptr[i], typeptr, soa.size());
        } else {
          std::tie(sigmaSquareds, epsilon24s) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeptr[i], typeptr, soa.size());
        }
      }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsSum)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        SoAFloatPrecision shift6 = h_shift6;
        SoAFloatPrecision sigmaSquared = h_sigmaSquared;
        SoAFloatPrecision epsilon24 = h_epsilon24;
        if constexpr (!isHomogeneous and useMixing) {
          sigmaSquared = sigmaSquareds[j];
          epsilon24 = epsilon24s[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

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
        const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // newton 3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          numKernelCallsN3Sum += mask;
        }

        if (calculateGlobals) {
          const SoAFloatPrecision virialx = drx * fx;
          const SoAFloatPrecision virialy = dry * fy;
          const SoAFloatPrecision virialz = drz * fz;
          const SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          SoAFloatPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                                           (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;

          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;

          if constexpr (countFLOPs) {
            numGlobalCalcsSum += mask;
          }
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsSum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) return;
    if (newton3) {
      dispatchSoAPair<true>(soa1, soa2);
    } else {
      dispatchSoAPair<false>(soa1, soa2);
    }
  }



  template <bool newton3>
  void dispatchSoAPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if constexpr (useMixing) {
      const auto *const typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
      const auto *const typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();
      if (_mixingTable.numTypes <= 1) {
        const auto t1 = typeptr1[0];
        const auto t2 = typeptr2[0];
        SoAFunctorPairImpl<true, newton3>(soa1, soa2,
                                          _mixingTable.template get<SigmaSquared>(t1, t2),
                                          _mixingTable.template get<Epsilon24>(t1, t2),
                                          applyShift ? _mixingTable.template get<Shift6>(t1, t2) : 0.);
      } else if (isCellHomogeneous(soa1) and isCellHomogeneous(soa2)) {
        const auto t1 = typeptr1[0];
        const auto t2 = typeptr2[0];
        SoAFunctorPairImpl<true, newton3>(soa1, soa2,
                                          _mixingTable.template get<SigmaSquared>(t1, t2),
                                          _mixingTable.template get<Epsilon24>(t1, t2),
                                          applyShift ? _mixingTable.template get<Shift6>(t1, t2) : 0.);
      } else {
        SoAFunctorPairImpl<false, newton3>(soa1, soa2, 0., 0., 0.);
      }
    } else {
      SoAFunctorPairImpl<true, newton3>(soa1, soa2, _sigmaSquared, _epsilon24, _shift6);
    }
  }

  template <bool isHomogeneous, bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                          SoAFloatPrecision h_sigmaSquared,
                          SoAFloatPrecision h_epsilon24,
                          SoAFloatPrecision h_shift6) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    [[maybe_unused]] PreloadCacheType preloadCache;
    if constexpr (!isHomogeneous and useMixing) {
      preloadCache = PreloadCacheType{_mixingTable.numTypes, soa2.size()};
    }

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision *sigmaSquareds = nullptr;
      const SoAFloatPrecision *epsilon24s = nullptr;
      const SoAFloatPrecision *shift6s = nullptr;
      if constexpr (!isHomogeneous and useMixing) {
        if constexpr (applyShift) {
          std::tie(sigmaSquareds, epsilon24s, shift6s) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeptr1[i], typeptr2, soa2.size());
        } else {
          std::tie(sigmaSquareds, epsilon24s) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeptr1[i], typeptr2, soa2.size());
        }
      }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        SoAFloatPrecision shift6 = h_shift6;
        SoAFloatPrecision sigmaSquared = h_sigmaSquared;
        SoAFloatPrecision epsilon24 = h_epsilon24;
        if constexpr (!isHomogeneous and useMixing) {
          sigmaSquared = sigmaSquareds[j];
          epsilon24 = epsilon24s[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

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
        const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

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

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          if constexpr (newton3) {
            numKernelCallsN3Sum += mask;
          } else {
            numKernelCallsNoN3Sum += mask;
          }
        }

        if constexpr (calculateGlobals) {
          SoAFloatPrecision virialx = drx * fx;
          SoAFloatPrecision virialy = dry * fy;
          SoAFloatPrecision virialz = drz * fz;
          SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          const SoAFloatPrecision energyFactor =
              (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
              (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;
          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;

          if constexpr (countFLOPs) {
            if constexpr (newton3) {
              numGlobalCalcsN3Sum += mask;
            } else {
              numGlobalCalcsNoN3Sum += mask;
            }
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

 public:
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        std::span<const size_t> neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if constexpr (useMixing) {
      if (_mixingTable.numTypes <= 1) {
        if (newton3) {
          SoAFunctorVerletImpl<true, true>(soa, indexFirst, neighborList);
        } else {
          SoAFunctorVerletImpl<false, true>(soa, indexFirst, neighborList);
        }
        return;
      }
    }
    if (newton3) {
      SoAFunctorVerletImpl<true, false>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false, false>(soa, indexFirst, neighborList);
    }
  }

  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    if (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0.;
    }
  }

  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 9>{Particle_T::AttributeNames::id,
                                                              Particle_T::AttributeNames::posX,
                                                              Particle_T::AttributeNames::posY,
                                                              Particle_T::AttributeNames::posZ,
                                                              Particle_T::AttributeNames::forceX,
                                                              Particle_T::AttributeNames::forceY,
                                                              Particle_T::AttributeNames::forceZ,
                                                              Particle_T::AttributeNames::typeId,
                                                              Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  constexpr static bool getMixing() { return useMixing; }

  void initTraversal() final {
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    if constexpr (calculateGlobals) {
      for (auto &data : _aosThreadDataGlobals) {
        data.setZero();
      }
    }
    if constexpr (countFLOPs) {
      for (auto &data : _aosThreadDataFLOPs) {
        data.setZero();
      }
    }
  }

  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum += data.virialSum;
      }
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;
      _potentialEnergySum /= 6.;
      _postProcessed = true;
    }
  }

  double getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get potential energy, because endTraversal was not called.");
    }
    return _potentialEnergySum;
  }

  double getVirial() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

  [[nodiscard]] size_t getNumFLOPs() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });
      const size_t numGlobalCalcsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsN3; });
      const size_t numGlobalCalcsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsNoN3; });

      constexpr size_t numFLOPsPerDistanceCall = 8;
      constexpr size_t numFLOPsPerN3KernelCall = 18;
      constexpr size_t numFLOPsPerNoN3KernelCall = 15;
      constexpr size_t numFLOPsPerN3GlobalCalc = applyShift ? 13 : 12;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = applyShift ? 9 : 8;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
    } else {
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numDistCallsAcc));
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

 private:
  template <bool newton3, bool isHomogeneous = false>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            std::span<const size_t> neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa.template begin<Particle_T::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision epsilon24 = _epsilon24;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;
    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

#ifdef __AVX512F__
    constexpr size_t vecsize = 16;
#else
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    const size_t typeI = (useMixing and not isHomogeneous) ? typeptr1[indexFirst] : 0;
    const SoAFloatPrecision *const __restrict__ sigmaSqRow =
        (useMixing and not isHomogeneous) ? _mixingTable.template row<LJMixParam::SigmaSquared>(typeI) : nullptr;
    const SoAFloatPrecision *const __restrict__ eps24Row =
        (useMixing and not isHomogeneous) ? _mixingTable.template row<LJMixParam::Epsilon24>(typeI) : nullptr;
    const SoAFloatPrecision *const __restrict__ shift6Row =
        (useMixing and not isHomogeneous and applyShift) ? _mixingTable.template row<LJMixParam::Shift6>(typeI) : nullptr;

    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
      }
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> sigmaSquareds;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> epsilon24s;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> shift6s;

#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          const size_t neighborIdx = neighborListPtr[joff + tmpj];
          xArr[tmpj] = xptr[neighborIdx];
          yArr[tmpj] = yptr[neighborIdx];
          zArr[tmpj] = zptr[neighborIdx];
          ownedStateArr[tmpj] = ownedStatePtr[neighborIdx];
          if constexpr (useMixing and not isHomogeneous) {
            const size_t typeJ = typeptr2[neighborIdx];
            sigmaSquareds[tmpj] = sigmaSqRow[typeJ];
            epsilon24s[tmpj] = eps24Row[typeJ];
            if constexpr (applyShift) {
              shift6s[tmpj] = shift6Row[typeJ];
            }
          }
        }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum) safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {
          if constexpr (useMixing and not isHomogeneous) {
            sigmaSquared = sigmaSquareds[j];
            epsilon24 = epsilon24s[j];
            if constexpr (applyShift) {
              shift6 = shift6s[j];
            }
          }

          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

          const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

          const SoAFloatPrecision invdr2 = 1. / dr2;
          const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
          const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
          const SoAFloatPrecision lj12 = lj6 * lj6;
          const SoAFloatPrecision lj12m6 = lj12 - lj6;
          const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

          const SoAFloatPrecision fx = drx * fac;
          const SoAFloatPrecision fy = dry * fac;
          const SoAFloatPrecision fz = drz * fac;

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if (newton3) {
            fxArr[j] = fx;
            fyArr[j] = fy;
            fzArr[j] = fz;
          }

          if constexpr (countFLOPs) {
            numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
            if constexpr (newton3) {
              numKernelCallsN3Sum += mask;
            } else {
              numKernelCallsNoN3Sum += mask;
            }
          }

          if (calculateGlobals) {
            SoAFloatPrecision virialx = drx * fx;
            SoAFloatPrecision virialy = dry * fy;
            SoAFloatPrecision virialz = drz * fz;
            SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

            const SoAFloatPrecision energyFactor =
                (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
            potentialEnergySum += potentialEnergy6 * energyFactor;
            virialSumX += virialx * energyFactor;
            virialSumY += virialy * energyFactor;
            virialSumZ += virialz * energyFactor;

            if constexpr (countFLOPs) {
              if constexpr (newton3) {
                numGlobalCalcsN3Sum += mask;
              } else {
                numGlobalCalcsNoN3Sum += mask;
              }
            }
          }
        }
        if (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborListPtr[joff + tmpj];
            fxptr[j] -= fxArr[tmpj];
            fyptr[j] -= fyArr[tmpj];
            fzptr[j] -= fzArr[tmpj];
          }
        }
      }
    }
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;
      if constexpr (useMixing and not isHomogeneous) {
        const size_t typeJ = typeptr2[j];
        sigmaSquared = sigmaSqRow[typeJ];
        epsilon24 = eps24Row[typeJ];
        if constexpr (applyShift) {
          shift6 = shift6Row[typeJ];
        }
      }

      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];

      const SoAFloatPrecision drx2 = drx * drx;
      const SoAFloatPrecision dry2 = dry * dry;
      const SoAFloatPrecision drz2 = drz * drz;

      const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

      if constexpr (countFLOPs) {
        numDistanceCalculationSum += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      const SoAFloatPrecision invdr2 = 1. / dr2;
      const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
      const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
      const SoAFloatPrecision lj12 = lj6 * lj6;
      const SoAFloatPrecision lj12m6 = lj12 - lj6;
      const SoAFloatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;

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

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          numKernelCallsN3Sum += 1;
        } else {
          numKernelCallsNoN3Sum += 1;
        }
      }

      if (calculateGlobals) {
        SoAFloatPrecision virialx = drx * fx;
        SoAFloatPrecision virialy = dry * fy;
        SoAFloatPrecision virialz = drz * fz;
        SoAFloatPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        const SoAFloatPrecision energyFactor =
            (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
            (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
        potentialEnergySum += potentialEnergy6 * energyFactor;
        virialSumX += virialx * energyFactor;
        virialSumY += virialy * energyFactor;
        virialSumZ += virialz * energyFactor;

        if constexpr (countFLOPs) {
          if constexpr (newton3) {
            ++numGlobalCalcsN3Sum;
          } else {
            ++numGlobalCalcsNoN3Sum;
          }
        }
      }
    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }

    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };

  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __remainingTo64{} {}

    void setZero() {
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numGlobalCalcsNoN3 = 0;
      numGlobalCalcsN3 = 0;
    }

    size_t numKernelCallsNoN3 = 0;
    size_t numKernelCallsN3 = 0;
    size_t numDistCalls = 0;
    size_t numGlobalCalcsN3 = 0;
    size_t numGlobalCalcsNoN3 = 0;

   private:
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;
  double _epsilon24{}, _sigmaSquared{}, _shift6 = 0;

  CompactMixingTable<2, 3, SoAFloatPrecision> _mixingTable{};

  double _potentialEnergySum;
  std::array<double, 3> _virialSum;

  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  bool _postProcessed;
};

}  // namespace mdLib
