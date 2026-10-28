/**
 * @file LJFunctorHWYSwissKnife.h
 * @brief Multi-regime Lennard-Jones Highway SIMD Functor adaptively choosing:
 *        1. Homogeneous register broadcast (0 memory loads) when cell(s) share a single type
 *        2. Row-hoist scalar gather for small heterogeneous cells (N <= 32)
 *        3. SoAMixingPreloadCache for large heterogeneous cells (N > 32)
 */

#pragma once

#include <hwy/highway.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "MixingTable.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/optRef.h"

namespace mdLib {

template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorHWYSwissKnife
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorHWYSwissKnife<Particle_T, applyShift, useMixing, useNewton3,
                                                                        calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using VectorizationPattern = autopas::VectorizationPatternOption::Value;

  static constexpr highway::ScalableTag<double> tag_double{};
  static constexpr highway::ScalableTag<int64_t> tag_long{};
  using VectorDouble = decltype(highway::Zero(tag_double));
  using VectorLong = decltype(highway::Zero(tag_long));
  using MaskDouble = decltype(highway::FirstN(tag_double, 1));

  enum LJMixParam : size_t { SigmaSquared = 0, Epsilon24 = 1, Shift6 = 2 };
  using CacheType = SoAMixingPreloadCache<double, applyShift ? 3 : 2, 4, 256>;

 public:
  LJFunctorHWYSwissKnife() = delete;

 private:
  explicit LJFunctorHWYSwissKnife(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorHWYSwissKnife<Particle_T, applyShift, useMixing, useNewton3,
                                                                   calculateGlobals, countFLOPs, relevantForTuning>>(
            cutoff),
        _cutoffSquare{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  explicit LJFunctorHWYSwissKnife(double cutoff) : LJFunctorHWYSwissKnife(cutoff, nullptr) {
    static_assert(not useMixing, "Mixing without ParticlePropertiesLibrary is not permitted.");
  }

  explicit LJFunctorHWYSwissKnife(double cutoff, ParticlePropertiesLibrary<double, size_t> &ppl)
      : LJFunctorHWYSwissKnife(cutoff, nullptr) {
    static_assert(useMixing, "Cannot use ParticlePropertiesLibrary with mixing disabled.");
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _mixingTable.build(nTypes, [&](size_t i, size_t j) -> std::array<double, 3> {
      const auto d = ppl.getLJMixingData(i, j);
      return {d.sigmaSquared, d.epsilon24, d.shift6};
    });
    if (nTypes > 0) {
      _sigmaSquared = _mixingTable.template get<SigmaSquared>(0, 0);
      _epsilon24 = _mixingTable.template get<Epsilon24>(0, 0);
      if constexpr (applyShift) _shift6 = _mixingTable.template get<Shift6>(0, 0);
    }
  }

  explicit LJFunctorHWYSwissKnife(double cutoff, std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> pplOpt)
      : LJFunctorHWYSwissKnife(cutoff, nullptr) {
    if constexpr (useMixing) {
      if (not pplOpt.has_value()) {
        autopas::utils::ExceptionHandler::exception("Mixing is enabled but no ParticlePropertiesLibrary was provided!");
      }
      auto &ppl = pplOpt.value().get();
      const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
      _mixingTable.build(nTypes, [&](size_t i, size_t j) -> std::array<double, 3> {
        const auto d = ppl.getLJMixingData(i, j);
        return {d.sigmaSquared, d.epsilon24, d.shift6};
      });
      if (nTypes > 0) {
        _sigmaSquared = _mixingTable.template get<SigmaSquared>(0, 0);
        _epsilon24 = _mixingTable.template get<Epsilon24>(0, 0);
        if constexpr (applyShift) _shift6 = _mixingTable.template get<Shift6>(0, 0);
      }
    } else {
      if (pplOpt.has_value()) {
        autopas::utils::ExceptionHandler::exception("Mixing is disabled but a ParticlePropertiesLibrary was provided!");
      }
    }
  }

  std::string getName() final { return "LJFunctorHWYSwissKnife"; }
  bool isRelevantForTuning() final { return relevantForTuning; }
  bool allowsNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Off; }
  bool allowsNonNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Only; }

  bool isVecPatternAllowed(const VectorizationPattern vecPattern) final {
    return vecPattern == VectorizationPattern::p1xVec;
  }

  void setVecPattern(const VectorizationPattern /*vecPattern*/) {}

  void setParticleProperties(double epsilon24, double sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    _mixingTable.resize(1);
    _mixingTable.channels[SigmaSquared].sbo[0] = sigmaSquared;
    _mixingTable.channels[Epsilon24].sbo[0] = epsilon24;
    _mixingTable.rebindPointers();
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

    const auto dr = i.getR() - j.getR();
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
    if (dr2 > _cutoffSquare) return;

    double sigSq = _sigmaSquared;
    double eps24 = _epsilon24;
    double sh6 = _shift6;
    if constexpr (useMixing) {
      sigSq = _mixingTable.template get<SigmaSquared>(i.getTypeId(), j.getTypeId());
      eps24 = _mixingTable.template get<Epsilon24>(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) sh6 = _mixingTable.template get<Shift6>(i.getTypeId(), j.getTypeId());
    }

    const double invdr2 = 1.0 / dr2;
    const double lj2 = sigSq * invdr2;
    const double lj6 = lj2 * lj2 * lj2;
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double fac = eps24 * (lj12 + lj12m6) * invdr2;
    const auto f = dr * fac;

    i.addF(f);
    if (newton3) j.subF(f);

    if constexpr (calculateGlobals) {
      const auto virial = dr * f;
      const double pot6 = eps24 * lj12m6 + sh6;
      const int tid = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[tid].potentialEnergySum += pot6;
        _aosThreadData[tid].virialSum += virial;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadData[tid].potentialEnergySum += pot6;
        _aosThreadData[tid].virialSum += virial;
      }
    }
  }

  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    if (soa.size() == 0) return;

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeIDptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    const size_t n = soa.size();
    const size_t vecLen = highway::Lanes(tag_double);

    // Check Regime 1: Homogeneous cell
    bool isHomogeneous = true;
    if constexpr (useMixing) {
      const size_t firstType = typeIDptr[0];
      for (size_t k = 1; k < n; ++k) {
        if (typeIDptr[k] != firstType) {
          isHomogeneous = false;
          break;
        }
      }
    }

    if (not useMixing or isHomogeneous) {
      // Regime 1: Register broadcast (zero memory loads for parameters)
      const double sigSq = useMixing ? _mixingTable.template get<SigmaSquared>(typeIDptr[0], typeIDptr[0]) : _sigmaSquared;
      const double eps24 = useMixing ? _mixingTable.template get<Epsilon24>(typeIDptr[0], typeIDptr[0]) : _epsilon24;
      double sh6 = _shift6;
      if constexpr (useMixing and applyShift) {
        sh6 = _mixingTable.template get<Shift6>(typeIDptr[0], typeIDptr[0]);
      }

      const auto sigmaSquareds = highway::Set(tag_double, sigSq);
      const auto epsilon24s = highway::Set(tag_double, eps24);
      const auto shift6s = highway::Set(tag_double, sh6);

      for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
        const auto x1 = highway::Set(tag_double, xPtr[i]);
        const auto y1 = highway::Set(tag_double, yPtr[i]);
        const auto z1 = highway::Set(tag_double, zPtr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = 0;
        const size_t limit = static_cast<size_t>(i);
        const size_t numFullSteps = limit / vecLen;
        const size_t fullLimit = numFullSteps * vecLen;

        for (; j < fullLimit; j += vecLen) {
          SoAKernelConst<true, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                      fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                      sigmaSquareds, epsilon24s, shift6s);
        }
        const size_t rest = limit - j;
        if (rest > 0) {
          SoAKernelConst<true, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                     fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                     sigmaSquareds, epsilon24s, shift6s);
        }

        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    } else if (n <= 8) {
      // Regime 2a: Small-N Direct Scalar Shortcut (bypasses vector remainder overhead and STLF)
      for (size_t i = 1; i < n; ++i) {
        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) continue;
        const double xi = xPtr[i], yi = yPtr[i], zi = zPtr[i];
        const size_t ti = typeIDptr[i];
        const double *const sigRow = _mixingTable.template row<SigmaSquared>(ti);
        const double *const epsRow = _mixingTable.template row<Epsilon24>(ti);
        const double *const shRow = applyShift ? _mixingTable.template row<Shift6>(ti) : nullptr;
        double fxi = 0., fyi = 0., fzi = 0.;
        for (size_t j = 0; j < i; ++j) {
          if (ownedStatePtr[j] == autopas::OwnershipState::dummy) continue;
          const double dx = xi - xPtr[j];
          const double dy = yi - yPtr[j];
          const double dz = zi - zPtr[j];
          const double dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 <= _cutoffSquare) {
            const size_t tj = typeIDptr[j];
            const double sigSq = sigRow[tj];
            const double eps24 = epsRow[tj];
            const double invdr2 = 1.0 / dr2;
            const double lj2 = sigSq * invdr2;
            const double lj4 = lj2 * lj2;
            const double lj6 = lj4 * lj2;
            const double lj12 = lj6 * lj6;
            const double lj12m6 = lj12 - lj6;
            const double fac = eps24 * (lj12 + lj12m6) * invdr2;
            const double fx = dx * fac;
            const double fy = dy * fac;
            const double fz = dz * fac;
            fxi += fx; fyi += fy; fzi += fz;
            fxPtr[j] -= fx; fyPtr[j] -= fy; fzPtr[j] -= fz;
            if constexpr (calculateGlobals) {
              const double pot6 = eps24 * lj12m6 + (applyShift ? shRow[tj] : 0.0);
              const int tid = autopas::autopas_get_thread_num();
              if (ownedStatePtr[i] == autopas::OwnershipState::owned) {
                _aosThreadData[tid].potentialEnergySum += pot6;
                _aosThreadData[tid].virialSum[0] += fx * dx;
                _aosThreadData[tid].virialSum[1] += fy * dy;
                _aosThreadData[tid].virialSum[2] += fz * dz;
              }
              if (ownedStatePtr[j] == autopas::OwnershipState::owned) {
                _aosThreadData[tid].potentialEnergySum += pot6;
                _aosThreadData[tid].virialSum[0] += fx * dx;
                _aosThreadData[tid].virialSum[1] += fy * dy;
                _aosThreadData[tid].virialSum[2] += fz * dz;
              }
            }
          }
        }
        fxPtr[i] += fxi; fyPtr[i] += fyi; fzPtr[i] += fzi;
      }
    } else if (n <= 16) {
      // Regime 2b: Heterogeneous Small-N (Row Hoisting gather)
      for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
        const auto x1 = highway::Set(tag_double, xPtr[i]);
        const auto y1 = highway::Set(tag_double, yPtr[i]);
        const auto z1 = highway::Set(tag_double, zPtr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        const auto *const rowSig = _mixingTable.template row<SigmaSquared>(typeIDptr[i]);
        const auto *const rowEps = _mixingTable.template row<Epsilon24>(typeIDptr[i]);
        const double *rowShift = nullptr;
        if constexpr (applyShift) {
          rowShift = _mixingTable.template row<Shift6>(typeIDptr[i]);
        }

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = 0;
        const size_t limit = static_cast<size_t>(i);
        const size_t numFullSteps = limit / vecLen;
        const size_t fullLimit = numFullSteps * vecLen;

        for (; j < fullLimit; j += vecLen) {
          SoAKernelGather<true, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, typeIDptr,
                                       fxPtr, fyPtr, fzPtr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                       virialSumZ, uPotSum, 0, rowSig, rowEps, rowShift);
        }
        const size_t rest = limit - j;
        if (rest > 0) {
          SoAKernelGather<true, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, typeIDptr,
                                      fxPtr, fyPtr, fzPtr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                      virialSumZ, uPotSum, rest, rowSig, rowEps, rowShift);
        }

        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    } else {
      // Regime 3: Heterogeneous Large-N (Preload Cache Streaming)
      CacheType mixingCache{_mixingTable.numTypes, n};

      for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
        const auto x1 = highway::Set(tag_double, xPtr[i]);
        const auto y1 = highway::Set(tag_double, yPtr[i]);
        const auto z1 = highway::Set(tag_double, zPtr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        const double *preloadedSig = nullptr;
        const double *preloadedEps = nullptr;
        const double *preloadedShift = nullptr;

        if constexpr (applyShift) {
          std::tie(preloadedSig, preloadedEps, preloadedShift) =
              mixingCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeIDptr[i], typeIDptr, n);
        } else {
          std::tie(preloadedSig, preloadedEps) =
              mixingCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeIDptr[i], typeIDptr, n);
        }

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = 0;
        const size_t limit = static_cast<size_t>(i);
        const size_t numFullSteps = limit / vecLen;
        const size_t fullLimit = numFullSteps * vecLen;

        for (; j < fullLimit; j += vecLen) {
          SoAKernelPreload<true, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                       fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                       preloadedSig, preloadedEps, preloadedShift);
        }
        const size_t rest = limit - j;
        if (rest > 0) {
          SoAKernelPreload<true, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                      fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                      preloadedSig, preloadedEps, preloadedShift);
        }

        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    }

    if constexpr (calculateGlobals) {
      accumulateGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) return;
    if (newton3) {
      SoAFunctorPairImpl<true, false>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false, false>(soa1, soa2);
    }
  }

  inline void SoAFunctorPairSorted(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   const autopas::SoASortingData &sortingData, bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) return;
    if (newton3) {
      SoAFunctorPairImpl<true, true>(soa1, soa2, sortingData);
    } else {
      SoAFunctorPairImpl<false, true>(soa1, soa2, sortingData);
    }
  }

  void initTraversal() final {
    _postProcessed = false;
    if constexpr (calculateGlobals) {
      for (auto &t : _aosThreadData) {
        t.potentialEnergySum = 0.;
        t.virialSum = {0., 0., 0.};
      }
    }
  }

  void endTraversal(bool newton3) final {
    if constexpr (calculateGlobals) {
      if (_postProcessed) return;
      _potentialEnergySum = 0.;
      _virialSum = {0., 0., 0.};
      for (const auto &t : _aosThreadData) {
        _potentialEnergySum += t.potentialEnergySum;
        _virialSum[0] += t.virialSum[0];
        _virialSum[1] += t.virialSum[1];
        _virialSum[2] += t.virialSum[2];
      }
      if (not newton3) {
        _potentialEnergySum *= 0.5;
        _virialSum[0] *= 0.5;
        _virialSum[1] *= 0.5;
        _virialSum[2] *= 0.5;
      }
      _postProcessed = true;
    }
  }

  double getPotentialEnergy() const { return _potentialEnergySum; }
  double getVirial() const { return _virialSum[0] + _virialSum[1] + _virialSum[2]; }

 private:
  // Core LJ interaction math given pre-resolved parameter vectors
  template <bool newton3, bool remainderJ>
  inline void computeInteraction(
      const size_t j, const MaskDouble &ownedMaskI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, const autopas::OwnershipState *const __restrict ownedState2Ptr,
      double *const __restrict fx2Ptr, double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restJ,
      const VectorDouble &sigmaSquareds,
      const VectorDouble &epsilon24s,
      const VectorDouble &shift6s) {

    VectorDouble x2, y2, z2;
    VectorLong ownedLong;
    if constexpr (remainderJ) {
      x2 = highway::LoadN(tag_double, &x2Ptr[j], restJ);
      y2 = highway::LoadN(tag_double, &y2Ptr[j], restJ);
      z2 = highway::LoadN(tag_double, &z2Ptr[j], restJ);
      ownedLong = highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedState2Ptr[j]), restJ);
    } else {
      x2 = highway::LoadU(tag_double, &x2Ptr[j]);
      y2 = highway::LoadU(tag_double, &y2Ptr[j]);
      z2 = highway::LoadU(tag_double, &z2Ptr[j]);
      ownedLong = highway::LoadU(tag_long, reinterpret_cast<const int64_t *>(&ownedState2Ptr[j]));
    }
    const MaskDouble ownedMaskJ = highway::RebindMask(tag_double, highway::Ne(ownedLong, highway::Zero(tag_long)));

    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);

    const auto drX2 = highway::Mul(drX, drX);
    const auto drY2 = highway::Mul(drY, drY);
    const auto drZ2 = highway::Mul(drZ, drZ);
    const auto dr2 = highway::Add(highway::Add(drX2, drY2), drZ2);

    const VectorDouble cutoffSquared = highway::Set(tag_double, _cutoffSquare);
    const auto dummyMask = highway::And(ownedMaskI, ownedMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, cutoffSquared);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    const auto invDr2 = highway::Div(highway::Set(tag_double, 1.0), dr2);
    const auto lj2 = highway::Mul(sigmaSquareds, invDr2);
    const auto lj4 = highway::Mul(lj2, lj2);
    const auto lj6 = highway::Mul(lj2, lj4);
    const auto lj12 = highway::Mul(lj6, lj6);
    const auto lj12m6 = highway::Sub(lj12, lj6);
    const auto lj12m6alj12 = highway::Add(lj12m6, lj12);
    const auto lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);
    const auto fac = highway::Mul(lj12m6alj12e, invDr2);

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    const VectorDouble fx = highway::Mul(drX, facMasked);
    const VectorDouble fy = highway::Mul(drY, facMasked);
    const VectorDouble fz = highway::Mul(drZ, facMasked);

    fxAcc = highway::Add(fxAcc, fx);
    fyAcc = highway::Add(fyAcc, fy);
    fzAcc = highway::Add(fzAcc, fz);

    if constexpr (calculateGlobals) {
      const auto virialX = highway::Mul(fx, drX);
      const auto virialY = highway::Mul(fy, drY);
      const auto virialZ = highway::Mul(fz, drZ);

      const auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      const auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::MaskedSet(tag_double, dummyMask, 1.0);
      if constexpr (newton3) {
        energyFactor = highway::Add(energyFactor, highway::MaskedSet(tag_double, dummyMask, 1.0));
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }

    if constexpr (newton3) {
      const VectorDouble fx2Old = remainderJ ? highway::LoadN(tag_double, &fx2Ptr[j], restJ)
                                             : highway::LoadU(tag_double, &fx2Ptr[j]);
      const VectorDouble fy2Old = remainderJ ? highway::LoadN(tag_double, &fy2Ptr[j], restJ)
                                             : highway::LoadU(tag_double, &fy2Ptr[j]);
      const VectorDouble fz2Old = remainderJ ? highway::LoadN(tag_double, &fz2Ptr[j], restJ)
                                             : highway::LoadU(tag_double, &fz2Ptr[j]);

      const VectorDouble fx2New = highway::Sub(fx2Old, fx);
      const VectorDouble fy2New = highway::Sub(fy2Old, fy);
      const VectorDouble fz2New = highway::Sub(fz2Old, fz);

      if constexpr (remainderJ) {
        highway::StoreN(fx2New, tag_double, &fx2Ptr[j], restJ);
        highway::StoreN(fy2New, tag_double, &fy2Ptr[j], restJ);
        highway::StoreN(fz2New, tag_double, &fz2Ptr[j], restJ);
      } else {
        highway::StoreU(fx2New, tag_double, &fx2Ptr[j]);
        highway::StoreU(fy2New, tag_double, &fy2Ptr[j]);
        highway::StoreU(fz2New, tag_double, &fz2Ptr[j]);
      }
    }
  }

  // Regime 1 Kernel: Constant pre-broadcast parameters
  template <bool newton3, bool remainderJ>
  inline void SoAKernelConst(
      const size_t j, const MaskDouble &ownedMaskI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, const autopas::OwnershipState *const __restrict ownedState2Ptr,
      double *const __restrict fx2Ptr, double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restJ,
      const VectorDouble &sigmaSquareds,
      const VectorDouble &epsilon24s,
      const VectorDouble &shift6s) {
    computeInteraction<newton3, remainderJ>(
        j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr, fx2Ptr, fy2Ptr, fz2Ptr,
        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restJ,
        sigmaSquareds, epsilon24s, shift6s);
  }

  // Regime 2 Kernel: Row-hoisting gather
  template <bool newton3, bool remainderJ>
  inline void SoAKernelGather(
      const size_t j, const MaskDouble &ownedMaskI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, const autopas::OwnershipState *const __restrict ownedState2Ptr,
      const size_t *const __restrict typeID2Ptr,
      double *const __restrict fx2Ptr, double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restJ,
      const double *const __restrict rowSig,
      const double *const __restrict rowEps,
      const double *const __restrict rowShift) {
    const size_t lanes = remainderJ ? restJ : highway::Lanes(tag_double);
    alignas(64) double sigBuf[8];
    alignas(64) double epsBuf[8];
    alignas(64) double shBuf[8];

    for (size_t k = 0; k < lanes; ++k) {
      const auto tJ = typeID2Ptr[j + k];
      sigBuf[k] = rowSig[tJ];
      epsBuf[k] = rowEps[tJ];
      if constexpr (applyShift) shBuf[k] = rowShift[tJ];
    }

    const auto sigmaSquareds = remainderJ ? highway::LoadN(tag_double, sigBuf, restJ) : highway::Load(tag_double, sigBuf);
    const auto epsilon24s = remainderJ ? highway::LoadN(tag_double, epsBuf, restJ) : highway::Load(tag_double, epsBuf);
    const auto shift6s = applyShift ? (remainderJ ? highway::LoadN(tag_double, shBuf, restJ) : highway::Load(tag_double, shBuf))
                                    : highway::Zero(tag_double);

    computeInteraction<newton3, remainderJ>(
        j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr, fx2Ptr, fy2Ptr, fz2Ptr,
        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restJ,
        sigmaSquareds, epsilon24s, shift6s);
  }

  // Regime 3 Kernel: Preload cache streaming
  template <bool newton3, bool remainderJ>
  inline void SoAKernelPreload(
      const size_t j, const MaskDouble &ownedMaskI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, const autopas::OwnershipState *const __restrict ownedState2Ptr,
      double *const __restrict fx2Ptr, double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restJ,
      const double *const __restrict preloadedSig,
      const double *const __restrict preloadedEps,
      const double *const __restrict preloadedShift) {
    VectorDouble sigmaSquareds;
    VectorDouble epsilon24s;
    VectorDouble shift6s;

    if constexpr (remainderJ) {
      sigmaSquareds = highway::LoadN(tag_double, &preloadedSig[j], restJ);
      epsilon24s = highway::LoadN(tag_double, &preloadedEps[j], restJ);
      if constexpr (applyShift) {
        shift6s = highway::LoadN(tag_double, &preloadedShift[j], restJ);
      } else {
        shift6s = highway::Zero(tag_double);
      }
    } else {
      sigmaSquareds = highway::LoadU(tag_double, &preloadedSig[j]);
      epsilon24s = highway::LoadU(tag_double, &preloadedEps[j]);
      if constexpr (applyShift) {
        shift6s = highway::LoadU(tag_double, &preloadedShift[j]);
      } else {
        shift6s = highway::Zero(tag_double);
      }
    }

    computeInteraction<newton3, remainderJ>(
        j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr, fx2Ptr, fy2Ptr, fz2Ptr,
        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restJ,
        sigmaSquareds, epsilon24s, shift6s);
  }

  template <bool newton3, bool sorted>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                 autopas::utils::optRef<const autopas::SoASortingData> sortingData = std::nullopt) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();

    const auto *const __restrict x1Ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeID1Ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();

    const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeID2Ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx1Ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1Ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1Ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();

    auto *const __restrict fx2Ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2Ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2Ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();

    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    const size_t vecLen = highway::Lanes(tag_double);

    // Check Regime 1: Homogeneous pair
    bool c1Homogeneous = true;
    bool c2Homogeneous = true;
    if constexpr (useMixing) {
      const size_t tA = typeID1Ptr[0];
      for (size_t k = 1; k < n1; ++k) {
        if (typeID1Ptr[k] != tA) { c1Homogeneous = false; break; }
      }
      const size_t tB = typeID2Ptr[0];
      for (size_t k = 1; k < n2; ++k) {
        if (typeID2Ptr[k] != tB) { c2Homogeneous = false; break; }
      }
    }

    if (not useMixing or (c1Homogeneous and c2Homogeneous)) {
      // Regime 1: Constant parameters for entire cell pair
      const double sigSq = useMixing ? _mixingTable.template get<SigmaSquared>(typeID1Ptr[0], typeID2Ptr[0]) : _sigmaSquared;
      const double eps24 = useMixing ? _mixingTable.template get<Epsilon24>(typeID1Ptr[0], typeID2Ptr[0]) : _epsilon24;
      double sh6 = _shift6;
      if constexpr (useMixing and applyShift) {
        sh6 = _mixingTable.template get<Shift6>(typeID1Ptr[0], typeID2Ptr[0]);
      }

      const auto sigmaSquareds = highway::Set(tag_double, sigSq);
      const auto epsilon24s = highway::Set(tag_double, eps24);
      const auto shift6s = highway::Set(tag_double, sh6);

      for (size_t i = 0; i < n1; ++i) {
        size_t jStart = 0;
        size_t jEnd = n2;
        if constexpr (sorted) {
          const auto &sd = sortingData->get();
          jEnd = sd.maxIndex[i];
          jStart = sd.minIndex[i];
          if (jStart >= jEnd) continue;
          jStart = jStart - (jStart % vecLen);
          if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) continue;
        }

        const auto x1 = highway::Set(tag_double, x1Ptr[i]);
        const auto y1 = highway::Set(tag_double, y1Ptr[i]);
        const auto z1 = highway::Set(tag_double, z1Ptr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr1[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = jStart;
        const size_t count = jEnd - jStart;
        const size_t fullSteps = count / vecLen;
        const size_t fullEnd = jStart + fullSteps * vecLen;

        if (newton3) {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelConst<true, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                        sigmaSquareds, epsilon24s, shift6s);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelConst<true, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                       fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                       sigmaSquareds, epsilon24s, shift6s);
          }
        } else {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelConst<false, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                         fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                         sigmaSquareds, epsilon24s, shift6s);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelConst<false, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                        sigmaSquareds, epsilon24s, shift6s);
          }
        }

        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    } else if (n2 <= 32) {
      // Regime 2: Heterogeneous Small-N2 (Row Hoisting gather)
      for (size_t i = 0; i < n1; ++i) {
        size_t jStart = 0;
        size_t jEnd = n2;
        if constexpr (sorted) {
          const auto &sd = sortingData->get();
          jEnd = sd.maxIndex[i];
          jStart = sd.minIndex[i];
          if (jStart >= jEnd) continue;
          jStart = jStart - (jStart % vecLen);
          if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) continue;
        }

        const auto x1 = highway::Set(tag_double, x1Ptr[i]);
        const auto y1 = highway::Set(tag_double, y1Ptr[i]);
        const auto z1 = highway::Set(tag_double, z1Ptr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr1[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        const auto *const rowSig = _mixingTable.template row<SigmaSquared>(typeID1Ptr[i]);
        const auto *const rowEps = _mixingTable.template row<Epsilon24>(typeID1Ptr[i]);
        const double *rowShift = nullptr;
        if constexpr (applyShift) {
          rowShift = _mixingTable.template row<Shift6>(typeID1Ptr[i]);
        }

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = jStart;
        const size_t count = jEnd - jStart;
        const size_t fullSteps = count / vecLen;
        const size_t fullEnd = jStart + fullSteps * vecLen;

        if (newton3) {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelGather<true, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, typeID2Ptr,
                                         fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                         virialSumZ, uPotSum, 0, rowSig, rowEps, rowShift);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelGather<true, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, typeID2Ptr,
                                        fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                        virialSumZ, uPotSum, rest, rowSig, rowEps, rowShift);
          }
        } else {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelGather<false, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, typeID2Ptr,
                                          fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                          virialSumZ, uPotSum, 0, rowSig, rowEps, rowShift);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelGather<false, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, typeID2Ptr,
                                         fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                         virialSumZ, uPotSum, rest, rowSig, rowEps, rowShift);
          }
        }

        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    } else {
      // Regime 3: Heterogeneous Large-N2 (Preload Cache Streaming)
      CacheType mixingCache{_mixingTable.numTypes, n2};

      for (size_t i = 0; i < n1; ++i) {
        size_t jStart = 0;
        size_t jEnd = n2;
        if constexpr (sorted) {
          const auto &sd = sortingData->get();
          jEnd = sd.maxIndex[i];
          jStart = sd.minIndex[i];
          if (jStart >= jEnd) continue;
          jStart = jStart - (jStart % vecLen);
          if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) continue;
        }

        const auto x1 = highway::Set(tag_double, x1Ptr[i]);
        const auto y1 = highway::Set(tag_double, y1Ptr[i]);
        const auto z1 = highway::Set(tag_double, z1Ptr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr1[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        const double *preloadedSig = nullptr;
        const double *preloadedEps = nullptr;
        const double *preloadedShift = nullptr;

        if constexpr (applyShift) {
          std::tie(preloadedSig, preloadedEps, preloadedShift) =
              mixingCache.template getPreloadedArrays<SigmaSquared, Epsilon24, Shift6>(
                  _mixingTable, typeID1Ptr[i], typeID2Ptr, n2);
        } else {
          std::tie(preloadedSig, preloadedEps) =
              mixingCache.template getPreloadedArrays<SigmaSquared, Epsilon24>(
                  _mixingTable, typeID1Ptr[i], typeID2Ptr, n2);
        }

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = jStart;
        const size_t count = jEnd - jStart;
        const size_t fullSteps = count / vecLen;
        const size_t fullEnd = jStart + fullSteps * vecLen;

        if (newton3) {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelPreload<true, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                         fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                         preloadedSig, preloadedEps, preloadedShift);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelPreload<true, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                        fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                        preloadedSig, preloadedEps, preloadedShift);
          }
        } else {
          for (; j < fullEnd; j += vecLen) {
            SoAKernelPreload<false, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                          fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0,
                                          preloadedSig, preloadedEps, preloadedShift);
          }
          const size_t rest = jEnd - j;
          if (rest > 0) {
            SoAKernelPreload<false, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                         fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest,
                                         preloadedSig, preloadedEps, preloadedShift);
          }
        }

        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
    }

    if constexpr (calculateGlobals) {
      accumulateGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  inline void accumulateGlobals(const VectorDouble &virialSumX, const VectorDouble &virialSumY,
                                const VectorDouble &virialSumZ, const VectorDouble &uPotSum) {
    const int tid = autopas::autopas_get_thread_num();
    _aosThreadData[tid].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
    _aosThreadData[tid].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
    _aosThreadData[tid].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    _aosThreadData[tid].potentialEnergySum += highway::ReduceSum(tag_double, uPotSum);
  }

  struct AoSThreadData {
    double potentialEnergySum{0.0};
    std::array<double, 3> virialSum{0.0, 0.0, 0.0};
  };

  double _cutoffSquare{0.};
  double _sigmaSquared{1.};
  double _epsilon24{24.};
  double _shift6{0.};
  double _potentialEnergySum{0.};
  std::array<double, 3> _virialSum{0., 0., 0.};
  bool _postProcessed{false};
  std::vector<AoSThreadData> _aosThreadData;

  CompactMixingTable<2, 3, double> _mixingTable;
};

}  // namespace mdLib
