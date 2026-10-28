/**
 * @file LJFunctorAVXNewMixing.h
 * @brief Simplified, zero-duplication LJFunctorAVX using CompactMixingTable + SoAMixingPreloadCache.
 */
#pragma once
#ifndef __AVX__
#pragma message "Requested to compile LJFunctorAVXNewMixing but AVX is not available!"
#else
#include <immintrin.h>
#endif

#include <array>
#include <span>
#include <tuple>

#include "MixingTable.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorAVXNewMixing
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorAVXNewMixing<Particle_T, applyShift, useMixing, useNewton3,
                                                                        calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  enum LJMixParam : size_t { SigmaSquared = 0, Epsilon24 = 1, Shift6 = 2 };
  using PreloadCacheType = SoAMixingPreloadCache<double, applyShift ? 3 : 2>;

 public:
  LJFunctorAVXNewMixing() = delete;

 private:
  explicit LJFunctorAVXNewMixing(double cutoff, void * /*dummy*/)
#ifdef __AVX__
      : autopas::PairwiseFunctor<Particle_T, LJFunctorAVXNewMixing<Particle_T, applyShift, useMixing, useNewton3,
                                                                   calculateGlobals, countFLOPs, relevantForTuning>>(
            cutoff),
        _cutoffSquared{_mm256_set1_pd(cutoff * cutoff)},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorAVXNewMixing with countFLOPs but FLOP counting is not implemented.");
    }
  }
#else
      : autopas::Functor<Particle_T, LJFunctorAVXNewMixing<Particle_T, applyShift, useMixing, useNewton3,
                                                           calculateGlobals, relevantForTuning>>(cutoff) {
    autopas::utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif
 public:
  explicit LJFunctorAVXNewMixing(double cutoff) : LJFunctorAVXNewMixing(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  explicit LJFunctorAVXNewMixing(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorAVXNewMixing(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    initMixingTableFromPPL(particlePropertiesLibrary);
  }

  void initMixingTableFromPPL(const ParticlePropertiesLibrary<double, size_t> &ppl) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _mixingTable.build(nTypes, [&](size_t i, size_t j) -> std::array<double, 3> {
      const auto data = ppl.getLJMixingData(i, j);
      return {data.sigmaSquared, data.epsilon24, data.shift6};
    });
  }

  std::string getName() final { return "LJFunctorAVXNewMixing"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  inline void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto sigmaSquared = _sigmaSquaredAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
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

    if (dr2 > _cutoffSquaredAoS) {
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
    if (calculateGlobals) {
      auto virial = dr * f;
      double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += virial;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += virial;
      }
    }
  }

  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final { SoAFunctorSingleImpl(soa); }

  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }



  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifdef __AVX__
    if (soa.size() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();

    [[maybe_unused]] PreloadCacheType preloadCache =
        useMixing ? PreloadCacheType{_mixingTable.numTypes, soa.size()} : PreloadCacheType{};

    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[i]));

      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);

      const double *sigmaSquaredsPtr = nullptr;
      const double *epsilon24sPtr = nullptr;
      const double *shift6sPtr = nullptr;
      if constexpr (useMixing) {
        if constexpr (applyShift) {
          std::tie(sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeIDptr[i], typeIDptr, soa.size());
        } else {
          std::tie(sigmaSquaredsPtr, epsilon24sPtr) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeIDptr[i], typeIDptr, soa.size());
        }
      }

      size_t j = 0;
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                               zptr, fxptr, fyptr, fzptr, sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr, fxacc, fyacc,
                               fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);
      }
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0) {
        SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                              zptr, fxptr, fyptr, fzptr, sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr, fxacc, fyacc,
                              fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);
      }

      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d hSumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
#endif
  }

  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
#ifdef __AVX__
    if (soa1.size() == 0 || soa2.size() == 0) return;

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

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();

    [[maybe_unused]] PreloadCacheType preloadCache =
        useMixing ? PreloadCacheType{_mixingTable.numTypes, soa2.size()} : PreloadCacheType{};

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr1[i]));

      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);
      const double *sigmaSquaredsPtr = nullptr;
      const double *epsilon24sPtr = nullptr;
      const double *shift6sPtr = nullptr;
      if constexpr (useMixing) {
        if constexpr (applyShift) {
          std::tie(sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeID1ptr[i], typeID2ptr, soa2.size());
        } else {
          std::tie(sigmaSquaredsPtr, epsilon24sPtr) =
              preloadCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeID1ptr[i], typeID2ptr, soa2.size());
        }
      }

      unsigned int j = 0;
      for (; j < (soa2.size() & ~(vecLength - 1)); j += 4) {
        SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                  y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr,
                                  fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);
      }
      const int rest = (int)(soa2.size() & (vecLength - 1));
      if (rest > 0)
        SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                 y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, sigmaSquaredsPtr, epsilon24sPtr, shift6sPtr,
                                 fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);

      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d hSumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
#endif
  }

#ifdef __AVX__
  template <bool newton3, bool remainderIsMasked>
  inline void SoAKernel(const size_t j, const __m256i ownedStateI, const int64_t *const __restrict ownedStatePtr2,
                        const __m256d &x1, const __m256d &y1, const __m256d &z1, const double *const __restrict x2ptr,
                        const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                        double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                        double *const __restrict fz2ptr, const double *const __restrict sigmaSquaredsPtr,
                        const double *const __restrict epsilon24sPtr, const double *const __restrict shift6sPtr,
                        __m256d &fxacc, __m256d &fyacc, __m256d &fzacc, __m256d *virialSumX, __m256d *virialSumY,
                        __m256d *virialSumZ, __m256d *potentialEnergySum, const unsigned int rest = 0) {
    __m256d epsilon24s = _epsilon24;
    __m256d sigmaSquareds = _sigmaSquared;
    __m256d shift6s = _shift6;
    if constexpr (useMixing) {
      epsilon24s = remainderIsMasked ? _mm256_maskload_pd(&epsilon24sPtr[j], _masks[rest - 1])
                                     : _mm256_load_pd(&epsilon24sPtr[j]);
      sigmaSquareds = remainderIsMasked ? _mm256_maskload_pd(&sigmaSquaredsPtr[j], _masks[rest - 1])
                                        : _mm256_load_pd(&sigmaSquaredsPtr[j]);
      if constexpr (applyShift) {
        shift6s = remainderIsMasked ? _mm256_maskload_pd(&shift6sPtr[j], _masks[rest - 1])
                                    : _mm256_load_pd(&shift6sPtr[j]);
      }
    }

    const __m256d x2 = remainderIsMasked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&x2ptr[j]);
    const __m256d y2 = remainderIsMasked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&y2ptr[j]);
    const __m256d z2 = remainderIsMasked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&z2ptr[j]);

    const __m256d drx = _mm256_sub_pd(x1, x2);
    const __m256d dry = _mm256_sub_pd(y1, y2);
    const __m256d drz = _mm256_sub_pd(z1, z2);

    const __m256d drx2 = _mm256_mul_pd(drx, drx);
    const __m256d dry2 = _mm256_mul_pd(dry, dry);
    const __m256d drz2 = _mm256_mul_pd(drz, drz);

    const __m256d dr2PART = _mm256_add_pd(drx2, dry2);
    const __m256d dr2 = _mm256_add_pd(dr2PART, drz2);

    const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffSquared, _CMP_LE_OS);

    const __m256i ownedStateJ = remainderIsMasked
                                    ? _mm256_castpd_si256(_mm256_maskload_pd(
                                          reinterpret_cast<double const *>(&ownedStatePtr2[j]), _masks[rest - 1]))
                                    : _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&ownedStatePtr2[j]));
    const __m256d dummyMask = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _zero, _CMP_NEQ_UQ);
    const __m256d cutoffDummyMask = _mm256_and_pd(cutoffMask, dummyMask);

    if (_mm256_movemask_pd(cutoffDummyMask) == 0) {
      return;
    }

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(sigmaSquareds, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, epsilon24s);
    const __m256d fac = _mm256_mul_pd(lj12m6alj12e, invdr2);

    const __m256d facMasked =
        remainderIsMasked ? _mm256_and_pd(fac, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
                          : _mm256_and_pd(fac, cutoffDummyMask);

    const __m256d fx = _mm256_mul_pd(drx, facMasked);
    const __m256d fy = _mm256_mul_pd(dry, facMasked);
    const __m256d fz = _mm256_mul_pd(drz, facMasked);

    fxacc = _mm256_add_pd(fxacc, fx);
    fyacc = _mm256_add_pd(fyacc, fy);
    fzacc = _mm256_add_pd(fzacc, fz);

    if constexpr (newton3) {
      const __m256d fx2 =
          remainderIsMasked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fx2ptr[j]);
      const __m256d fy2 =
          remainderIsMasked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fy2ptr[j]);
      const __m256d fz2 =
          remainderIsMasked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fz2ptr[j]);

      const __m256d fx2new = _mm256_sub_pd(fx2, fx);
      const __m256d fy2new = _mm256_sub_pd(fy2, fy);
      const __m256d fz2new = _mm256_sub_pd(fz2, fz);

      remainderIsMasked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new)
                        : _mm256_storeu_pd(&fx2ptr[j], fx2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new)
                        : _mm256_storeu_pd(&fy2ptr[j], fy2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new)
                        : _mm256_storeu_pd(&fz2ptr[j], fz2new);
    }

    if constexpr (calculateGlobals) {
      const __m256d virialX = _mm256_mul_pd(fx, drx);
      const __m256d virialY = _mm256_mul_pd(fy, dry);
      const __m256d virialZ = _mm256_mul_pd(fz, drz);

      const __m256d potentialEnergy = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const __m256d potentialEnergyMasked =
          remainderIsMasked
              ? _mm256_and_pd(potentialEnergy, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
              : _mm256_and_pd(potentialEnergy, cutoffDummyMask);

      __m256d ownedMaskI =
          _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateI), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
      __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskI);
      if constexpr (newton3) {
        __m256d ownedMaskJ =
            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
        energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskJ));
      }
      *potentialEnergySum = wrapperFMA(energyFactor, potentialEnergyMasked, *potentialEnergySum);
      *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
      *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
      *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
    }
  }
#endif

 public:
  inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                               std::span<const size_t> neighborList,
                               bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

 private:
  template <bool newton3>
  inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                   std::span<const size_t> neighborList) {
#ifdef __AVX__
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();
    __m256d fxacc = _mm256_setzero_pd();
    __m256d fyacc = _mm256_setzero_pd();
    __m256d fzacc = _mm256_setzero_pd();

    const auto x1 = _mm256_broadcast_sd(&xptr[indexFirst]);
    const auto y1 = _mm256_broadcast_sd(&yptr[indexFirst]);
    const auto z1 = _mm256_broadcast_sd(&zptr[indexFirst]);
    __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[indexFirst]));

    alignas(64) std::array<double, vecLength> x2tmp{};
    alignas(64) std::array<double, vecLength> y2tmp{};
    alignas(64) std::array<double, vecLength> z2tmp{};
    alignas(64) std::array<double, vecLength> fx2tmp{};
    alignas(64) std::array<double, vecLength> fy2tmp{};
    alignas(64) std::array<double, vecLength> fz2tmp{};
    alignas(64) std::array<double, vecLength> sigmaSq2tmp{};
    alignas(64) std::array<double, vecLength> eps242tmp{};
    alignas(64) std::array<double, vecLength> shift62tmp{};
    alignas(64) std::array<autopas::OwnershipState, vecLength> ownedStates2tmp{};

    const size_t typeI = useMixing ? typeIDptr[indexFirst] : 0;
    const double *const __restrict__ sigmaSqRow =
        useMixing ? _mixingTable.template row<LJMixParam::SigmaSquared>(typeI) : nullptr;
    const double *const __restrict__ eps24Row =
        useMixing ? _mixingTable.template row<LJMixParam::Epsilon24>(typeI) : nullptr;
    const double *const __restrict__ shift6Row =
        (useMixing and applyShift) ? _mixingTable.template row<LJMixParam::Shift6>(typeI) : nullptr;

    size_t j = 0;
    for (; j < (neighborList.size() & ~(vecLength - 1)); j += vecLength) {
      for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
        const size_t neighborIdx = neighborList[j + vecIndex];
        x2tmp[vecIndex] = xptr[neighborIdx];
        y2tmp[vecIndex] = yptr[neighborIdx];
        z2tmp[vecIndex] = zptr[neighborIdx];
        if constexpr (newton3) {
          fx2tmp[vecIndex] = fxptr[neighborIdx];
          fy2tmp[vecIndex] = fyptr[neighborIdx];
          fz2tmp[vecIndex] = fzptr[neighborIdx];
        }
        if constexpr (useMixing) {
          const size_t typeJ = typeIDptr[neighborIdx];
          sigmaSq2tmp[vecIndex] = sigmaSqRow[typeJ];
          eps242tmp[vecIndex] = eps24Row[typeJ];
          if constexpr (applyShift) {
            shift62tmp[vecIndex] = shift6Row[typeJ];
          }
        }
        ownedStates2tmp[vecIndex] = ownedStatePtr[neighborIdx];
      }

      SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                                x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                                sigmaSq2tmp.data(), eps242tmp.data(), shift62tmp.data(), fxacc, fyacc, fzacc,
                                &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
          fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
          fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
          fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
        }
      }
    }

    const auto rest = static_cast<int>(neighborList.size() & (vecLength - 1));
    if (rest > 0) {
      for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
        const size_t neighborIdx = neighborList[j + vecIndex];
        x2tmp[vecIndex] = xptr[neighborIdx];
        y2tmp[vecIndex] = yptr[neighborIdx];
        z2tmp[vecIndex] = zptr[neighborIdx];
        if constexpr (newton3) {
          fx2tmp[vecIndex] = fxptr[neighborIdx];
          fy2tmp[vecIndex] = fyptr[neighborIdx];
          fz2tmp[vecIndex] = fzptr[neighborIdx];
        }
        if constexpr (useMixing) {
          const size_t typeJ = typeIDptr[neighborIdx];
          sigmaSq2tmp[vecIndex] = sigmaSqRow[typeJ];
          eps242tmp[vecIndex] = eps24Row[typeJ];
          if constexpr (applyShift) {
            shift62tmp[vecIndex] = shift6Row[typeJ];
          }
        }
        ownedStates2tmp[vecIndex] = ownedStatePtr[neighborIdx];
      }

      SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                               x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                               sigmaSq2tmp.data(), eps242tmp.data(), shift62tmp.data(), fxacc, fyacc, fzacc,
                               &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
          fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
          fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
          fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
        }
      }
    }

    const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
    const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

    const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
    const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

    const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
    const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

    const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
    const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

    const double sumfx = sumfxfyVEC[0];
    const double sumfy = sumfxfyVEC[1];
    const double sumfz = _mm_cvtsd_f64(sumfzVEC);

    fxptr[indexFirst] += sumfx;
    fyptr[indexFirst] += sumfy;
    fzptr[indexFirst] += sumfz;

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d hSumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
#endif
  }

 public:
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
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
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

  void setParticleProperties(double epsilon24, double sigmaSquared) {
#ifdef __AVX__
    _epsilon24 = _mm256_set1_pd(epsilon24);
    _sigmaSquared = _mm256_set1_pd(sigmaSquared);
    if constexpr (applyShift) {
      _shift6 = _mm256_set1_pd(
          ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared[0]));
    } else {
      _shift6 = _mm256_setzero_pd();
    }
#endif

    _epsilon24AoS = epsilon24;
    _sigmaSquaredAoS = sigmaSquared;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquaredAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:
#ifdef __AVX__
  inline __m256d wrapperFMA(const __m256d &factorA, const __m256d &factorB, const __m256d &summandC) {
#ifdef __FMA__
    return _mm256_fmadd_pd(factorA, factorB, summandC);
#elif __AVX__
    const __m256d tmp = _mm256_mul_pd(factorA, factorB);
    return _mm256_add_pd(summandC, tmp);
#else
    return __m256d();
#endif
  }
#endif

  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

#ifdef __AVX__
  const __m256d _zero{_mm256_set1_pd(0.)};
  const __m256d _one{_mm256_set1_pd(1.)};
  const __m256i _vindex = _mm256_set_epi64x(0, 1, 3, 4);
  const __m256i _masks[3]{
      _mm256_set_epi64x(0, 0, 0, -1),
      _mm256_set_epi64x(0, 0, -1, -1),
      _mm256_set_epi64x(0, -1, -1, -1),
  };
  const __m256i _ownedStateDummyMM256i{0x0};
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(autopas::OwnershipState::owned))};
  const __m256d _cutoffSquared{};
  __m256d _shift6 = _mm256_setzero_pd();
  __m256d _epsilon24{};
  __m256d _sigmaSquared{};
#endif

  const double _cutoffSquaredAoS = 0;
  double _epsilon24AoS{}, _sigmaSquaredAoS{}, _shift6AoS = 0;

  CompactMixingTable<2, 3, double> _mixingTable{};

  double _potentialEnergySum;
  std::array<double, 3> _virialSum;
  std::vector<AoSThreadData> _aosThreadData;
  bool _postProcessed;

  constexpr static size_t vecLength = 4;
};

}  // namespace mdLib
