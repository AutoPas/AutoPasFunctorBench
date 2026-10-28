/**
 * @file LJFunctorHWYOnTheFly1D.h
 * Functor variant that computes Lorentz-Berthelot mixing on-the-fly using Highway SIMD,
 * but projects per-type parameters (halfSigma, sqrtEpsilon) into temporary 1D contiguous
 * arrays at cell entry, rather than storing them permanently in the particle SoA.
 */

#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <optional>
#include <functional>

#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"

#include <hwy/highway.h>

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;

template <class Particle_T, bool applyShift = false, bool useMixing = true,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorHWYOnTheFly1D
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorHWYOnTheFly1D<Particle_T, applyShift, useMixing,
                                                                          useNewton3, calculateGlobals, countFLOPs,
                                                                          relevantForTuning>> {
  using VectorDouble = highway::Vec<highway::ScalableTag<double>>;
  using MaskDouble = highway::Mask<highway::ScalableTag<double>>;
  using VectorLong = highway::Vec<highway::ScalableTag<int64_t>>;
  using MaskLong = highway::Mask<highway::ScalableTag<int64_t>>;

  static constexpr highway::ScalableTag<double> tag_double{};
  static constexpr highway::ScalableTag<int64_t> tag_long{};

  struct AoSThreadData {
    double potentialEnergySum{0.};
    std::array<double, 3> virialSum{0., 0., 0.};
  };

 public:
  using SoAArraysType = typename Particle_T::SoAArraysType;

  explicit LJFunctorHWYOnTheFly1D(double cutoff)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorHWYOnTheFly1D<Particle_T, applyShift, useMixing,
                                                                    useNewton3, calculateGlobals, countFLOPs,
                                                                    relevantForTuning>>(cutoff),
        _cutoffSquare{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

  explicit LJFunctorHWYOnTheFly1D(double cutoff, ParticlePropertiesLibrary<double, size_t> &ppl)
      : LJFunctorHWYOnTheFly1D(cutoff) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _halfSigmaPerType.resize(nTypes);
    _sqrtEpsilonPerType.resize(nTypes);
    for (size_t t = 0; t < nTypes; ++t) {
      _halfSigmaPerType[t] = ppl.getSigma(t) * 0.5;
      _sqrtEpsilonPerType[t] = std::sqrt(ppl.getEpsilon(t));
    }
    if (nTypes > 0) {
      _sigmaSquared = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
      _epsilon24 = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
    }
  }

  explicit LJFunctorHWYOnTheFly1D(double cutoff, std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> pplOpt)
      : LJFunctorHWYOnTheFly1D(cutoff) {
    if (pplOpt.has_value()) {
      auto &ppl = pplOpt.value().get();
      const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
      _halfSigmaPerType.resize(nTypes);
      _sqrtEpsilonPerType.resize(nTypes);
      for (size_t t = 0; t < nTypes; ++t) {
        _halfSigmaPerType[t] = ppl.getSigma(t) * 0.5;
        _sqrtEpsilonPerType[t] = std::sqrt(ppl.getEpsilon(t));
      }
      if (nTypes > 0) {
        _sigmaSquared = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
        _epsilon24 = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
      }
    }
  }

  std::string getName() final { return "LJFunctorHWYOnTheFly1D"; }
  bool isRelevantForTuning() final { return relevantForTuning; }
  bool allowsNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Off; }
  bool allowsNonNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Only; }

  double getPotentialEnergy() const { return _potentialEnergySum; }
  double getVirial() const { return _virialSum[0] + _virialSum[1] + _virialSum[2]; }

  bool isVecPatternAllowed(const VectorizationPattern vecPattern) final {
    return vecPattern == VectorizationPattern::p1xVec;
  }
  void setVecPattern(const VectorizationPattern) {}

  void setParticleProperties(double epsilon24, double sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    _halfSigmaPerType = {std::sqrt(sigmaSquared) * 0.5};
    _sqrtEpsilonPerType = {std::sqrt(epsilon24 / 24.0)};
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

    const size_t ti = i.getTypeId();
    const size_t tj = j.getTypeId();
    const double sig = _halfSigmaPerType[ti] + _halfSigmaPerType[tj];
    const double sigSq = sig * sig;
    const double eps24 = 24.0 * _sqrtEpsilonPerType[ti] * _sqrtEpsilonPerType[tj];

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
      const double pot6 = eps24 * lj12m6;
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
    const size_t n = soa.size();
    if (n == 0) return;

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeIDptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    // Small-N Scalar Shortcut: Bypasses vector remainder overhead and horizontal reductions
    if (n <= 8) {
      for (size_t i = 1; i < n; ++i) {
        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) continue;
        const double xi = xPtr[i], yi = yPtr[i], zi = zPtr[i];
        const size_t ti = typeIDptr[i];
        const double halfSigI = _halfSigmaPerType[ti];
        const double sqrtEpsI = _sqrtEpsilonPerType[ti];
        double fxi = 0., fyi = 0., fzi = 0.;
        for (size_t j = 0; j < i; ++j) {
          if (ownedStatePtr[j] == autopas::OwnershipState::dummy) continue;
          const double dx = xi - xPtr[j];
          const double dy = yi - yPtr[j];
          const double dz = zi - zPtr[j];
          const double dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 <= _cutoffSquare) {
            const size_t tj = typeIDptr[j];
            const double sig = halfSigI + _halfSigmaPerType[tj];
            const double sigSq = sig * sig;
            const double eps24 = 24.0 * sqrtEpsI * _sqrtEpsilonPerType[tj];
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
              const double pot6 = eps24 * lj12m6;
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
      return;
    }

    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    const size_t vecLen = highway::Lanes(tag_double);

    // Fast-path for Homogeneous single type (zero parameter loads)
    if (_halfSigmaPerType.size() <= 1) {
      const auto sigSq = highway::Set(tag_double, _sigmaSquared);
      const auto eps24 = highway::Set(tag_double, _epsilon24);

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
          computeKernel<true, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                     fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0, sigSq, eps24);
        }
        const size_t rest = limit - j;
        if (rest > 0) {
          computeKernel<true, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                    fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest, sigSq, eps24);
        }

        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
      return;
    }

    // Heterogeneous Separable 1D On-The-Fly Mixing: Unpack parameters into contiguous stack buffer
    alignas(64) double halfSigmaStack[256];
    alignas(64) double sqrtEpsStack[256];
    double *halfSigmaPtr = halfSigmaStack;
    double *sqrtEpsPtr = sqrtEpsStack;
    std::vector<double, autopas::AlignedAllocator<double, 64>> heapHalfSig, heapSqrtEps;
    if (n > 256) {
      heapHalfSig.resize(n);
      heapSqrtEps.resize(n);
      halfSigmaPtr = heapHalfSig.data();
      sqrtEpsPtr = heapSqrtEps.data();
    }
    for (size_t k = 0; k < n; ++k) {
      const auto t = typeIDptr[k];
      halfSigmaPtr[k] = _halfSigmaPerType[t];
      sqrtEpsPtr[k] = _sqrtEpsilonPerType[t];
    }

    const auto vec24 = highway::Set(tag_double, 24.0);

    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
      const auto x1 = highway::Set(tag_double, xPtr[i]);
      const auto y1 = highway::Set(tag_double, yPtr[i]);
      const auto z1 = highway::Set(tag_double, zPtr[i]);
      const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
      const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

      const auto halfSigI = highway::Set(tag_double, halfSigmaPtr[i]);
      const auto sqrtEpsI = highway::Set(tag_double, sqrtEpsPtr[i]);

      VectorDouble fxAcc = highway::Zero(tag_double);
      VectorDouble fyAcc = highway::Zero(tag_double);
      VectorDouble fzAcc = highway::Zero(tag_double);

      size_t j = 0;
      const size_t limit = static_cast<size_t>(i);
      const size_t numFullSteps = limit / vecLen;
      const size_t fullLimit = numFullSteps * vecLen;

      for (; j < fullLimit; j += vecLen) {
        const auto halfSigJ = highway::LoadU(tag_double, &halfSigmaPtr[j]);
        const auto sqrtEpsJ = highway::LoadU(tag_double, &sqrtEpsPtr[j]);
        const auto sig = highway::Add(halfSigI, halfSigJ);
        const auto sigSq = highway::Mul(sig, sig);
        const auto eps24 = highway::Mul(vec24, highway::Mul(sqrtEpsI, sqrtEpsJ));

        computeKernel<true, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                   fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0, sigSq, eps24);
      }
      const size_t rest = limit - j;
      if (rest > 0) {
        const auto halfSigJ = highway::LoadN(tag_double, &halfSigmaPtr[j], rest);
        const auto sqrtEpsJ = highway::LoadN(tag_double, &sqrtEpsPtr[j], rest);
        const auto sig = highway::Add(halfSigI, halfSigJ);
        const auto sigSq = highway::Mul(sig, sig);
        const auto eps24 = highway::Mul(vec24, highway::Mul(sqrtEpsI, sqrtEpsJ));

        computeKernel<true, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                  fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest, sigSq, eps24);
      }

      fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
      fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
      fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
    }
  }

  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) final {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();
    if (n1 == 0 or n2 == 0) return;

    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> &soa1, autopas::SoAView<SoAArraysType> &soa2) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();

    const auto *const __restrict x1Ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedState1Ptr = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeID1ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx1Ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1Ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1Ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedState2Ptr = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx2Ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2Ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2Ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();

    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    const size_t vecLen = highway::Lanes(tag_double);

    // Fast-path for Homogeneous single type (zero parameter loads)
    if (_halfSigmaPerType.size() <= 1) {
      const auto sigSq = highway::Set(tag_double, _sigmaSquared);
      const auto eps24 = highway::Set(tag_double, _epsilon24);
      const size_t numFullSteps = n2 / vecLen;
      const size_t fullLimit = numFullSteps * vecLen;
      const size_t rest = n2 - fullLimit;

      for (size_t i = 0; i < n1; ++i) {
        const auto x1 = highway::Set(tag_double, x1Ptr[i]);
        const auto y1 = highway::Set(tag_double, y1Ptr[i]);
        const auto z1 = highway::Set(tag_double, z1Ptr[i]);
        const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedState1Ptr[i]));
        const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

        VectorDouble fxAcc = highway::Zero(tag_double);
        VectorDouble fyAcc = highway::Zero(tag_double);
        VectorDouble fzAcc = highway::Zero(tag_double);

        size_t j = 0;
        for (; j < fullLimit; j += vecLen) {
          computeKernel<newton3, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr,
                                       fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                       virialSumZ, uPotSum, 0, sigSq, eps24);
        }
        if (rest > 0) {
          computeKernel<newton3, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr,
                                      fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                      virialSumZ, uPotSum, rest, sigSq, eps24);
        }

        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
      }
      return;
    }

    // Unpack 1D parameter projections on-the-fly for both cells
    alignas(64) double halfSigmaStack1[256], sqrtEpsStack1[256];
    alignas(64) double halfSigmaStack2[256], sqrtEpsStack2[256];
    double *halfSig1Ptr = halfSigmaStack1, *sqrtEps1Ptr = sqrtEpsStack1;
    double *halfSig2Ptr = halfSigmaStack2, *sqrtEps2Ptr = sqrtEpsStack2;
    std::vector<double, autopas::AlignedAllocator<double, 64>> hSig1, hEps1, hSig2, hEps2;
    if (n1 > 256) {
      hSig1.resize(n1); hEps1.resize(n1);
      halfSig1Ptr = hSig1.data(); sqrtEps1Ptr = hEps1.data();
    }
    for (size_t k = 0; k < n1; ++k) {
      const auto t = typeID1ptr[k];
      halfSig1Ptr[k] = _halfSigmaPerType[t];
      sqrtEps1Ptr[k] = _sqrtEpsilonPerType[t];
    }
    if (n2 > 256) {
      hSig2.resize(n2); hEps2.resize(n2);
      halfSig2Ptr = hSig2.data(); sqrtEps2Ptr = hEps2.data();
    }
    for (size_t k = 0; k < n2; ++k) {
      const auto t = typeID2ptr[k];
      halfSig2Ptr[k] = _halfSigmaPerType[t];
      sqrtEps2Ptr[k] = _sqrtEpsilonPerType[t];
    }

    const auto vec24 = highway::Set(tag_double, 24.0);
    const size_t numFullSteps = n2 / vecLen;
    const size_t fullLimit = numFullSteps * vecLen;
    const size_t rest = n2 - fullLimit;

    for (size_t i = 0; i < n1; ++i) {
      const auto x1 = highway::Set(tag_double, x1Ptr[i]);
      const auto y1 = highway::Set(tag_double, y1Ptr[i]);
      const auto z1 = highway::Set(tag_double, z1Ptr[i]);
      const auto ownedStateLongI = highway::Set(tag_long, static_cast<int64_t>(ownedState1Ptr[i]));
      const MaskDouble ownedMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateLongI, highway::Zero(tag_long)));

      const auto halfSigI = highway::Set(tag_double, halfSig1Ptr[i]);
      const auto sqrtEpsI = highway::Set(tag_double, sqrtEps1Ptr[i]);

      VectorDouble fxAcc = highway::Zero(tag_double);
      VectorDouble fyAcc = highway::Zero(tag_double);
      VectorDouble fzAcc = highway::Zero(tag_double);

      size_t j = 0;
      for (; j < fullLimit; j += vecLen) {
        const auto halfSigJ = highway::LoadU(tag_double, &halfSig2Ptr[j]);
        const auto sqrtEpsJ = highway::LoadU(tag_double, &sqrtEps2Ptr[j]);
        const auto sig = highway::Add(halfSigI, halfSigJ);
        const auto sigSq = highway::Mul(sig, sig);
        const auto eps24 = highway::Mul(vec24, highway::Mul(sqrtEpsI, sqrtEpsJ));

        computeKernel<newton3, false>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr,
                                     fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                     virialSumZ, uPotSum, 0, sigSq, eps24);
      }
      if (rest > 0) {
        const auto halfSigJ = highway::LoadN(tag_double, &halfSig2Ptr[j], rest);
        const auto sqrtEpsJ = highway::LoadN(tag_double, &sqrtEps2Ptr[j], rest);
        const auto sig = highway::Add(halfSigI, halfSigJ);
        const auto sigSq = highway::Mul(sig, sig);
        const auto eps24 = highway::Mul(vec24, highway::Mul(sqrtEpsI, sqrtEpsJ));

        computeKernel<newton3, true>(j, ownedMaskI, x1, y1, z1, x2Ptr, y2Ptr, z2Ptr, ownedState2Ptr,
                                    fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                    virialSumZ, uPotSum, rest, sigSq, eps24);
      }

      fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
      fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
      fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
    }
  }

  template <bool newton3, bool remainderJ>
  [[gnu::always_inline]] inline void computeKernel(
      const size_t j, const MaskDouble &ownedMaskI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, const autopas::OwnershipState *const __restrict ownedState2Ptr,
      double *const __restrict fx2Ptr, double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restJ, const VectorDouble &sigmaSquareds, const VectorDouble &epsilon24s) {
    const VectorDouble x2 = remainderJ ? highway::LoadN(tag_double, &x2Ptr[j], restJ)
                                       : highway::LoadU(tag_double, &x2Ptr[j]);
    const VectorDouble y2 = remainderJ ? highway::LoadN(tag_double, &y2Ptr[j], restJ)
                                       : highway::LoadU(tag_double, &y2Ptr[j]);
    const VectorDouble z2 = remainderJ ? highway::LoadN(tag_double, &z2Ptr[j], restJ)
                                       : highway::LoadU(tag_double, &z2Ptr[j]);

    VectorLong ownedLong;
    if constexpr (remainderJ) {
      ownedLong = highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedState2Ptr[j]), restJ);
    } else {
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

    if (highway::AllFalse(tag_double, cutoffDummyMask)) return;

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

  const double _cutoffSquare;
  double _potentialEnergySum;
  std::array<double, 3> _virialSum;
  bool _postProcessed;
  std::vector<AoSThreadData> _aosThreadData{};

  std::vector<double> _halfSigmaPerType{};
  std::vector<double> _sqrtEpsilonPerType{};

  double _sigmaSquared{0.};
  double _epsilon24{0.};
};

}  // namespace mdLib
