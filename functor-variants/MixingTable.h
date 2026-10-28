/**
 * @file MixingTable.h
 * @brief Generic, cache-line-aligned Structure-of-Arrays (SoA) mixing table and stack preload cache.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"

namespace mdLib {

/**
 * Generic N-body Structure-of-Arrays (SoA) mixing table with Small-Buffer Optimization (SBO).
 *
 * @tparam Arity Interaction order: 2 for pairwise (M x M), 3 for triwise (M x M x M).
 * @tparam NumParams Number of mixed floating-point quantities stored per type combination.
 * @tparam FloatType Floating-point precision (default: double).
 * @tparam MaxStaticTypes Maximum number of particle types stored inline via SBO before heap fallback (default: 8).
 */
template <size_t Arity, size_t NumParams, typename FloatType = double, size_t MaxStaticTypes = 8>
class CompactMixingTable {
  static_assert(Arity == 2 or Arity == 3, "CompactMixingTable currently supports Arity == 2 or Arity == 3.");

 public:
  static constexpr size_t MaxTypes = MaxStaticTypes;
  static constexpr size_t StaticCapacity =
      (Arity == 2) ? (MaxStaticTypes * MaxStaticTypes) : (MaxStaticTypes * MaxStaticTypes * MaxStaticTypes);

  struct Channel {
    alignas(64) std::array<FloatType, StaticCapacity> sbo{};
    std::vector<FloatType, autopas::AlignedAllocator<FloatType, 64>> dyn{};
    const FloatType *__restrict__ ptr{nullptr};
  };

  size_t numTypes{0};
  std::array<Channel, NumParams> channels{};

  CompactMixingTable() noexcept { rebindPointers(); }

  CompactMixingTable(const CompactMixingTable &other)
      : numTypes{other.numTypes}, channels{other.channels} {
    rebindPointers();
  }

  CompactMixingTable(CompactMixingTable &&other) noexcept
      : numTypes{other.numTypes}, channels{std::move(other.channels)} {
    rebindPointers();
  }

  CompactMixingTable &operator=(const CompactMixingTable &other) {
    if (this != &other) {
      numTypes = other.numTypes;
      channels = other.channels;
      rebindPointers();
    }
    return *this;
  }

  CompactMixingTable &operator=(CompactMixingTable &&other) noexcept {
    if (this != &other) {
      numTypes = other.numTypes;
      channels = std::move(other.channels);
      rebindPointers();
    }
    return *this;
  }

  void rebindPointers() noexcept {
    const bool useSBO = (numTypes <= MaxStaticTypes);
    for (size_t p = 0; p < NumParams; ++p) {
      channels[p].ptr = useSBO ? channels[p].sbo.data() : channels[p].dyn.data();
    }
  }

  void resize(size_t nTypes) {
    numTypes = nTypes;
    if (nTypes > MaxStaticTypes) {
      const size_t totalEntries = (Arity == 2) ? (nTypes * nTypes) : (nTypes * nTypes * nTypes);
      for (size_t p = 0; p < NumParams; ++p) {
        channels[p].dyn.resize(totalEntries);
      }
    }
    rebindPointers();
  }

  template <typename RuleFn>
  void build(size_t nTypes, RuleFn &&mixingRule) {
    resize(nTypes);
    if constexpr (Arity == 2) {
      for (size_t i = 0; i < nTypes; ++i) {
        for (size_t j = 0; j < nTypes; ++j) {
          const std::array<FloatType, NumParams> vals = mixingRule(i, j);
          const size_t idx = i * nTypes + j;
          for (size_t p = 0; p < NumParams; ++p) {
            writeEntry(p, idx, vals[p]);
          }
        }
      }
    } else if constexpr (Arity == 3) {
      for (size_t i = 0; i < nTypes; ++i) {
        for (size_t j = 0; j < nTypes; ++j) {
          for (size_t k = 0; k < nTypes; ++k) {
            const std::array<FloatType, NumParams> vals = mixingRule(i, j, k);
            const size_t idx = (i * nTypes + j) * nTypes + k;
            for (size_t p = 0; p < NumParams; ++p) {
              writeEntry(p, idx, vals[p]);
            }
          }
        }
      }
    }
  }

  // Direct 2D element getter (for AoS)
  template <size_t P>
  [[gnu::always_inline]] FloatType get(size_t i, size_t j) const noexcept {
    static_assert(Arity == 2 and P < NumParams);
    return channels[P].ptr[i * numTypes + j];
  }

  // Direct 3D element getter (for 3-body AoS)
  template <size_t P>
  [[gnu::always_inline]] FloatType get(size_t i, size_t j, size_t k) const noexcept {
    static_assert(Arity == 3 and P < NumParams);
    return channels[P].ptr[(i * numTypes + j) * numTypes + k];
  }

  // 1D row pointer for pairwise functors
  template <size_t P>
  [[gnu::always_inline]] const FloatType *row(size_t i) const noexcept {
    static_assert(Arity == 2 and P < NumParams);
    return channels[P].ptr + i * numTypes;
  }

  // 1D row pointer for triwise functors
  template <size_t P>
  [[gnu::always_inline]] const FloatType *row(size_t i, size_t j) const noexcept {
    static_assert(Arity == 3 and P < NumParams);
    return channels[P].ptr + (i * numTypes + j) * numTypes;
  }

 private:
  void writeEntry(size_t p, size_t flatIdx, FloatType val) noexcept {
    if (numTypes <= MaxStaticTypes) {
      channels[p].sbo[flatIdx] = val;
    } else {
      channels[p].dyn[flatIdx] = val;
    }
  }
};

/**
 * Stack-allocated SBO cache that lazily preloads 64B-aligned contiguous parameter arrays
 * once per unique outer particle type (O(M * N) instead of O(N^2)).
 *
 * Works identically for auto-vectorized (#pragma omp simd) and manually vectorized (AVX2 / Highway) functors.
 */
template <typename FloatType, size_t NumActiveParams, size_t MaxCachedTypes = 4, size_t MaxStackParticles = 256>
class SoAMixingPreloadCache {
 public:
  struct alignas(64) AlignedBuffer {
    FloatType data[MaxStackParticles];
  };

  SoAMixingPreloadCache() noexcept = default;

  SoAMixingPreloadCache(size_t numTypes, size_t numParticlesJ) noexcept
      : _useStack{numTypes <= MaxCachedTypes and numParticlesJ <= MaxStackParticles} {
    if (not _useStack) {
      _heapStorage = std::make_unique<HeapStorage>();
      _heapStorage->buffers.resize(NumActiveParams);
      for (size_t p = 0; p < NumActiveParams; ++p) {
        _heapStorage->buffers[p].resize(numTypes);
        for (size_t t = 0; t < numTypes; ++t) {
          _heapStorage->buffers[p][t].resize(numParticlesJ);
        }
      }
      _heapStorage->populated.assign(numTypes, false);
    }
  }

  /**
   * Returns a std::tuple of 64-byte-aligned contiguous pointers [param0_arr, param1_arr, ...] for outer type `typeI`.
   * Lazily populates the row on the first occurrence of `typeI` in the cell.
   * Can be directly unpacked via structured binding: auto [p0, p1] = ...
   * or assigned to existing pointers via std::tie(p0, p1) = ...
   */
  template <size_t... ParamIndices, typename TableType>
  [[gnu::always_inline]] auto getPreloadedArrays(
      const TableType &table, size_t typeI, const size_t *const __restrict__ typeptrJ, size_t numParticlesJ) noexcept {
    static_assert(sizeof...(ParamIndices) == NumActiveParams);

    if (_useStack) {
      if ((_stackPopulatedMask & (1u << typeI)) == 0) {
        _stackPopulatedMask |= (1u << typeI);
        populateRow<0, ParamIndices...>(table, typeI, typeptrJ, numParticlesJ);
      }
      return [&]<size_t... Is>(std::index_sequence<Is...>) {
        return std::make_tuple(
            static_cast<const FloatType *__restrict__>(_stackBuf[Is][typeI].data)...);
      }(std::make_index_sequence<NumActiveParams>{});
    } else {
      if (not _heapStorage->populated[typeI]) {
        _heapStorage->populated[typeI] = true;
        populateRow<0, ParamIndices...>(table, typeI, typeptrJ, numParticlesJ);
      }
      return [&]<size_t... Is>(std::index_sequence<Is...>) {
        return std::make_tuple(
            static_cast<const FloatType *__restrict__>(_heapStorage->buffers[Is][typeI].data())...);
      }(std::make_index_sequence<NumActiveParams>{});
    }
  }

 private:
  template <size_t SlotIdx, size_t FirstParam, size_t... RestParams, typename TableType>
  [[gnu::always_inline]] void populateRow(const TableType &table, size_t typeI,
                                          const size_t *const __restrict__ typeptrJ, size_t numParticlesJ) noexcept {
    FloatType *__restrict__ dst =
        _useStack ? _stackBuf[SlotIdx][typeI].data : _heapStorage->buffers[SlotIdx][typeI].data();
    const FloatType *const __restrict__ srcRow = table.template row<FirstParam>(typeI);

#pragma omp simd
    for (size_t j = 0; j < numParticlesJ; ++j) {
      dst[j] = srcRow[typeptrJ[j]];
    }

    if constexpr (sizeof...(RestParams) > 0) {
      populateRow<SlotIdx + 1, RestParams...>(table, typeI, typeptrJ, numParticlesJ);
    }
  }

  struct HeapStorage {
    std::vector<std::vector<std::vector<FloatType, autopas::AlignedAllocator<FloatType, 64>>>> buffers;
    std::vector<bool> populated;
  };

  bool _useStack{true};
  uint32_t _stackPopulatedMask{0};
  alignas(64) AlignedBuffer _stackBuf[NumActiveParams][MaxCachedTypes];
  std::unique_ptr<HeapStorage> _heapStorage{nullptr};
};

}  // namespace mdLib
