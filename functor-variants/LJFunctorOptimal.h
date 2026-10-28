/**
 * @file LJFunctorOptimal.h
 * @brief Auto-vectorized Lennard-Jones Functor utilizing the Optimal 1D On-The-Fly Mixing strategy.
 *
 * Key Architectural Highlights:
 * 1. Zero Particle Bloat: Particles only store standard attributes (positions, typeId, ownershipState).
 * 2. 1D Contiguous Parameter Projection: For separable mixing (Lorentz-Berthelot), per-type parameters
 *    are projected into 64-byte aligned stack buffers at cell entry (O(N) setup, ~5 ns).
 * 3. Multi-Regime Hardware Specialization:
 *    - Regime 1: Homogeneous Single-Type Fast Path (Scalar broadcast, zero parameter array reads).
 *    - Regime 2: Small-N Direct Scalar Shortcut (N <= 8 in SoASingle, avoids SIMD peeling/masking overhead).
 *    - Regime 3: Auto-vectorized 1D Contiguous Streaming (#pragma omp simd with contiguous array loads).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

template <class Particle_T, bool applyShift = false, bool useMixing = true,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorOptimal
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorOptimal<Particle_T, applyShift, useMixing, useNewton3,
                                                                   calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  struct AoSThreadData {
    double potentialEnergySum{0.};
    std::array<double, 3> virialSum{0., 0., 0.};
  };

 public:
  LJFunctorOptimal() = delete;

  explicit LJFunctorOptimal(double cutoff)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorOptimal<Particle_T, applyShift, useMixing, useNewton3,
                                                             calculateGlobals, countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{static_cast<SoAFloatPrecision>(cutoff * cutoff)},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

  explicit LJFunctorOptimal(double cutoff, ParticlePropertiesLibrary<double, size_t> &ppl)
      : LJFunctorOptimal(cutoff) {
    initFromPPL(ppl);
  }

  explicit LJFunctorOptimal(double cutoff,
                            std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> pplOpt)
      : LJFunctorOptimal(cutoff) {
    if (pplOpt.has_value()) {
      initFromPPL(pplOpt.value().get());
    }
  }

  std::string getName() final { return "LJFunctorOptimal"; }
  bool isRelevantForTuning() final { return relevantForTuning; }
  bool allowsNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Off; }
  bool allowsNonNewton3() final { return useNewton3 != autopas::FunctorN3Modes::Newton3Only; }

  double getPotentialEnergy() const { return _potentialEnergySum; }
  double getVirial() const { return _virialSum[0] + _virialSum[1] + _virialSum[2]; }

  void setParticleProperties(double epsilon24, double sigmaSquared) {
    _epsilon24 = static_cast<SoAFloatPrecision>(epsilon24);
    _sigmaSquared = static_cast<SoAFloatPrecision>(sigmaSquared);
    _halfSigmaPerType = {static_cast<SoAFloatPrecision>(std::sqrt(sigmaSquared) * 0.5)};
    _sqrtEpsilonPerType = {static_cast<SoAFloatPrecision>(std::sqrt(epsilon24 / 24.0))};
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

  // ==========================================
  // AoS Functor
  // ==========================================
  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) return;

    const auto dr = i.getR() - j.getR();
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
    if (dr2 > _cutoffSquared) return;

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

  // ==========================================
  // SoAFunctorSingle (Self-Interaction)
  // ==========================================
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    const size_t n = soa.size();
    if (n == 0) return;

    const auto *const __restrict x = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict owned = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    auto *const __restrict fx = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto cutoffSq = _cutoffSquared;

    // -------------------------------------------------------------------------
    // Regime 2: Small-N Direct Scalar Shortcut (N <= 8)
    // Avoids SIMD peel/tail loop branching on tiny loop counts
    // -------------------------------------------------------------------------
    if (n <= 8) {
      for (size_t i = 1; i < n; ++i) {
        if (owned[i] == autopas::OwnershipState::dummy) continue;
        const auto xi = x[i], yi = y[i], zi = z[i];
        const size_t ti = typeptr[i];
        const auto halfSigI = _halfSigmaPerType[ti];
        const auto sqrtEpsI = _sqrtEpsilonPerType[ti];
        SoAFloatPrecision fxi = 0., fyi = 0., fzi = 0.;

        for (size_t j = 0; j < i; ++j) {
          if (owned[j] == autopas::OwnershipState::dummy) continue;
          const auto dx = xi - x[j], dy = yi - y[j], dz = zi - z[j];
          const auto dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 <= cutoffSq) {
            const size_t tj = typeptr[j];
            const auto sig = halfSigI + _halfSigmaPerType[tj];
            const auto sigSq = sig * sig;
            const auto eps24 = static_cast<SoAFloatPrecision>(24.0) * sqrtEpsI * _sqrtEpsilonPerType[tj];
            const auto invdr2 = static_cast<SoAFloatPrecision>(1.0) / dr2;
            const auto lj2 = sigSq * invdr2;
            const auto lj4 = lj2 * lj2;
            const auto lj6 = lj4 * lj2;
            const auto lj12 = lj6 * lj6;
            const auto lj12m6 = lj12 - lj6;
            const auto fac = eps24 * (lj12 + lj12m6) * invdr2;
            const auto dfx = dx * fac, dfy = dy * fac, dfz = dz * fac;
            fxi += dfx; fyi += dfy; fzi += dfz;
            fx[j] -= dfx; fy[j] -= dfy; fz[j] -= dfz;

            if constexpr (calculateGlobals) {
              const auto pot6 = eps24 * lj12m6;
              const int tid = autopas::autopas_get_thread_num();
              if (owned[i] == autopas::OwnershipState::owned) {
                _aosThreadData[tid].potentialEnergySum += pot6;
                _aosThreadData[tid].virialSum[0] += dfx * dx;
                _aosThreadData[tid].virialSum[1] += dfy * dy;
                _aosThreadData[tid].virialSum[2] += dfz * dz;
              }
              if (owned[j] == autopas::OwnershipState::owned) {
                _aosThreadData[tid].potentialEnergySum += pot6;
                _aosThreadData[tid].virialSum[0] += dfx * dx;
                _aosThreadData[tid].virialSum[1] += dfy * dy;
                _aosThreadData[tid].virialSum[2] += dfz * dz;
              }
            }
          }
        }
        fx[i] += fxi; fy[i] += fyi; fz[i] += fzi;
      }
      return;
    }

    // -------------------------------------------------------------------------
    // Check cell homogeneity (All particles have the same type)
    // -------------------------------------------------------------------------
    bool isHomogeneous = true;
    const size_t firstType = typeptr[0];
    for (size_t k = 1; k < n; ++k) {
      if (typeptr[k] != firstType) { isHomogeneous = false; break; }
    }

    SoAFloatPrecision peSum = 0., virX = 0., virY = 0., virZ = 0.;

    // -------------------------------------------------------------------------
    // Regime 1: Homogeneous Fast Path (Scalar parameters, zero array loads)
    // -------------------------------------------------------------------------
    if (isHomogeneous) {
      const auto s2 = static_cast<SoAFloatPrecision>(4.0) * _halfSigmaPerType[firstType] * _halfSigmaPerType[firstType];
      const auto e24 = static_cast<SoAFloatPrecision>(24.0) * _sqrtEpsilonPerType[firstType] * _sqrtEpsilonPerType[firstType];

      for (unsigned int i = 0; i < n; ++i) {
        if (owned[i] == autopas::OwnershipState::dummy) continue;
        SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
        for (unsigned int j = i + 1; j < n; ++j) {
          const auto drx = x[i] - x[j], dry = y[i] - y[j], drz = z[i] - z[j];
          const auto dr2 = drx * drx + dry * dry + drz * drz;
          const bool mask = dr2 <= cutoffSq and owned[j] != autopas::OwnershipState::dummy;

          const auto invdr2 = static_cast<SoAFloatPrecision>(1.0) / dr2;
          const auto lj2 = s2 * invdr2;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : static_cast<SoAFloatPrecision>(0.0);

          const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
          fxacc += dfx; fyacc += dfy; fzacc += dfz;
          fx[j] -= dfx; fy[j] -= dfy; fz[j] -= dfz;

          if constexpr (calculateGlobals) {
            if (mask) {
              const auto pe = e24 * (lj12 - lj6);
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
        const auto tid = autopas::autopas_get_thread_num();
        _aosThreadData[tid].potentialEnergySum += peSum;
        _aosThreadData[tid].virialSum[0] += virX;
        _aosThreadData[tid].virialSum[1] += virY;
        _aosThreadData[tid].virialSum[2] += virZ;
      }
      return;
    }

    // -------------------------------------------------------------------------
    // Regime 3: Heterogeneous 1D On-The-Fly SIMD Streaming
    // Project parameters into contiguous stack arrays
    // -------------------------------------------------------------------------
    alignas(64) SoAFloatPrecision halfSigBuf[512];
    alignas(64) SoAFloatPrecision sqrtEpsBuf[512];
    SoAFloatPrecision *halfSig = halfSigBuf;
    SoAFloatPrecision *sqrtEps = sqrtEpsBuf;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision, 64>> heapHalfSig, heapSqrtEps;
    if (n > 512) {
      heapHalfSig.resize(n); heapSqrtEps.resize(n);
      halfSig = heapHalfSig.data(); sqrtEps = heapSqrtEps.data();
    }
    for (size_t k = 0; k < n; ++k) {
      const auto t = typeptr[k];
      halfSig[k] = _halfSigmaPerType[t];
      sqrtEps[k] = _sqrtEpsilonPerType[t];
    }

    for (unsigned int i = 0; i < n; ++i) {
      if (owned[i] == autopas::OwnershipState::dummy) continue;
      SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;
      const auto halfSigI = halfSig[i];
      const auto sqrtEpsI = sqrtEps[i];

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
      for (unsigned int j = i + 1; j < n; ++j) {
        const auto sig = halfSigI + halfSig[j];
        const auto s2 = sig * sig;
        const auto e24 = static_cast<SoAFloatPrecision>(24.0) * sqrtEpsI * sqrtEps[j];

        const auto drx = x[i] - x[j], dry = y[i] - y[j], drz = z[i] - z[j];
        const auto dr2 = drx * drx + dry * dry + drz * drz;
        const bool mask = dr2 <= cutoffSq and owned[j] != autopas::OwnershipState::dummy;

        const auto invdr2 = static_cast<SoAFloatPrecision>(1.0) / dr2;
        const auto lj2 = s2 * invdr2;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : static_cast<SoAFloatPrecision>(0.0);

        const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
        fxacc += dfx; fyacc += dfy; fzacc += dfz;
        fx[j] -= dfx; fy[j] -= dfy; fz[j] -= dfz;

        if constexpr (calculateGlobals) {
          if (mask) {
            const auto pe = e24 * (lj12 - lj6);
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
      const auto tid = autopas::autopas_get_thread_num();
      _aosThreadData[tid].potentialEnergySum += peSum;
      _aosThreadData[tid].virialSum[0] += virX;
      _aosThreadData[tid].virialSum[1] += virY;
      _aosThreadData[tid].virialSum[2] += virZ;
    }
  }

  // ==========================================
  // SoAFunctorPair (Cell Pair Interaction)
  // ==========================================
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) final {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();
    if (n1 == 0 or n2 == 0) return;

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

    const auto cutoffSq = _cutoffSquared;

    // Check pair homogeneity
    bool isPairHomogeneous = true;
    const size_t firstType = type1[0];
    for (size_t k = 0; k < n1; ++k) {
      if (type1[k] != firstType) { isPairHomogeneous = false; break; }
    }
    if (isPairHomogeneous) {
      for (size_t k = 0; k < n2; ++k) {
        if (type2[k] != firstType) { isPairHomogeneous = false; break; }
      }
    }

    SoAFloatPrecision peSum = 0., virX = 0., virY = 0., virZ = 0.;

    // Regime 1: Homogeneous Pair Fast Path
    if (isPairHomogeneous) {
      const auto s2 = static_cast<SoAFloatPrecision>(4.0) * _halfSigmaPerType[firstType] * _halfSigmaPerType[firstType];
      const auto e24 = static_cast<SoAFloatPrecision>(24.0) * _sqrtEpsilonPerType[firstType] * _sqrtEpsilonPerType[firstType];

      for (unsigned int i = 0; i < n1; ++i) {
        if (owned1[i] == autopas::OwnershipState::dummy) continue;
        SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
        for (unsigned int j = 0; j < n2; ++j) {
          const auto drx = x1[i] - x2[j], dry = y1[i] - y2[j], drz = z1[i] - z2[j];
          const auto dr2 = drx * drx + dry * dry + drz * drz;
          const bool mask = dr2 <= cutoffSq and owned2[j] != autopas::OwnershipState::dummy;

          const auto invdr2 = static_cast<SoAFloatPrecision>(1.0) / dr2;
          const auto lj2 = s2 * invdr2;
          const auto lj6 = lj2 * lj2 * lj2;
          const auto lj12 = lj6 * lj6;
          const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : static_cast<SoAFloatPrecision>(0.0);

          const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
          fxacc += dfx; fyacc += dfy; fzacc += dfz;
          if constexpr (useNewton3 != autopas::FunctorN3Modes::Newton3Off) {
            if (newton3) {
              fx2[j] -= dfx; fy2[j] -= dfy; fz2[j] -= dfz;
            }
          }

          if constexpr (calculateGlobals) {
            if (mask) {
              const auto pe = e24 * (lj12 - lj6);
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
        const auto tid = autopas::autopas_get_thread_num();
        _aosThreadData[tid].potentialEnergySum += peSum;
        _aosThreadData[tid].virialSum[0] += virX;
        _aosThreadData[tid].virialSum[1] += virY;
        _aosThreadData[tid].virialSum[2] += virZ;
      }
      return;
    }

    // Regime 3: Heterogeneous 1D On-The-Fly Mixing
    alignas(64) SoAFloatPrecision halfSig1Buf[512], sqrtEps1Buf[512];
    alignas(64) SoAFloatPrecision halfSig2Buf[512], sqrtEps2Buf[512];
    SoAFloatPrecision *halfSig1 = halfSig1Buf, *sqrtEps1 = sqrtEps1Buf;
    SoAFloatPrecision *halfSig2 = halfSig2Buf, *sqrtEps2 = sqrtEps2Buf;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision, 64>> hSig1, hEps1, hSig2, hEps2;
    if (n1 > 512) {
      hSig1.resize(n1); hEps1.resize(n1);
      halfSig1 = hSig1.data(); sqrtEps1 = hEps1.data();
    }
    for (size_t k = 0; k < n1; ++k) {
      const auto t = type1[k];
      halfSig1[k] = _halfSigmaPerType[t];
      sqrtEps1[k] = _sqrtEpsilonPerType[t];
    }
    if (n2 > 512) {
      hSig2.resize(n2); hEps2.resize(n2);
      halfSig2 = hSig2.data(); sqrtEps2 = hEps2.data();
    }
    for (size_t k = 0; k < n2; ++k) {
      const auto t = type2[k];
      halfSig2[k] = _halfSigmaPerType[t];
      sqrtEps2[k] = _sqrtEpsilonPerType[t];
    }

    for (unsigned int i = 0; i < n1; ++i) {
      if (owned1[i] == autopas::OwnershipState::dummy) continue;
      SoAFloatPrecision fxacc = 0., fyacc = 0., fzacc = 0.;
      const auto halfSigI = halfSig1[i];
      const auto sqrtEpsI = sqrtEps1[i];

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, peSum, virX, virY, virZ)
      for (unsigned int j = 0; j < n2; ++j) {
        const auto sig = halfSigI + halfSig2[j];
        const auto s2 = sig * sig;
        const auto e24 = static_cast<SoAFloatPrecision>(24.0) * sqrtEpsI * sqrtEps2[j];

        const auto drx = x1[i] - x2[j], dry = y1[i] - y2[j], drz = z1[i] - z2[j];
        const auto dr2 = drx * drx + dry * dry + drz * drz;
        const bool mask = dr2 <= cutoffSq and owned2[j] != autopas::OwnershipState::dummy;

        const auto invdr2 = static_cast<SoAFloatPrecision>(1.0) / dr2;
        const auto lj2 = s2 * invdr2;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto fac = mask ? e24 * (lj12 + lj12 - lj6) * invdr2 : static_cast<SoAFloatPrecision>(0.0);

        const auto dfx = drx * fac, dfy = dry * fac, dfz = drz * fac;
        fxacc += dfx; fyacc += dfy; fzacc += dfz;
        if constexpr (useNewton3 != autopas::FunctorN3Modes::Newton3Off) {
          if (newton3) {
            fx2[j] -= dfx; fy2[j] -= dfy; fz2[j] -= dfz;
          }
        }

        if constexpr (calculateGlobals) {
          if (mask) {
            const auto pe = e24 * (lj12 - lj6);
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
      const auto tid = autopas::autopas_get_thread_num();
      _aosThreadData[tid].potentialEnergySum += peSum;
      _aosThreadData[tid].virialSum[0] += virX;
      _aosThreadData[tid].virialSum[1] += virY;
      _aosThreadData[tid].virialSum[2] += virZ;
    }
  }

 private:
  void initFromPPL(ParticlePropertiesLibrary<double, size_t> &ppl) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _halfSigmaPerType.resize(nTypes);
    _sqrtEpsilonPerType.resize(nTypes);
    for (size_t t = 0; t < nTypes; ++t) {
      _halfSigmaPerType[t] = static_cast<SoAFloatPrecision>(ppl.getSigma(t) * 0.5);
      _sqrtEpsilonPerType[t] = static_cast<SoAFloatPrecision>(std::sqrt(ppl.getEpsilon(t)));
    }
    if (nTypes > 0) {
      _sigmaSquared = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
      _epsilon24 = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
    }
  }

  const SoAFloatPrecision _cutoffSquared;
  double _potentialEnergySum;
  std::array<double, 3> _virialSum;
  bool _postProcessed;
  std::vector<AoSThreadData> _aosThreadData{};

  std::vector<SoAFloatPrecision> _halfSigmaPerType{};
  std::vector<SoAFloatPrecision> _sqrtEpsilonPerType{};

  SoAFloatPrecision _sigmaSquared{0.};
  SoAFloatPrecision _epsilon24{0.};
};

}  // namespace mdLib
