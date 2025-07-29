#ifndef _TREETOOLS_CLUSTERTABLE_SIMD_H
#define _TREETOOLS_CLUSTERTABLE_SIMD_H

#include "ClusterTable.h"  // Include your existing header

// SIMD feature detection
#ifdef __AVX2__
#include <immintrin.h>
#define SIMD_WIDTH_32 8  // 8 x int32 per AVX2 vector
#define HAS_AVX2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define SIMD_WIDTH_32 4  // 4 x int32 per SSE2 vector
#define HAS_SSE2 1
#else
#define SIMD_WIDTH_32 1
#endif

namespace TreeTools {

// Simple SIMD utility functions
namespace simd_utils {

// Feature detection functions
inline bool has_avx2_support() {
#ifdef HAS_AVX2
  return true;
#else
  return false;
#endif
}

inline bool has_sse2_support() {
#ifdef HAS_SSE2
  return true;
#else
  return false;
#endif
}

// SIMD-optimized threshold counting for int32 arrays
inline int32 count_above_threshold(const int32* data, int32 size, int32 threshold) {
#ifdef HAS_AVX2
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m256i thresh_vec = _mm256_set1_epi32(threshold - 1);  // For >= comparison
  __m256i count_vec = _mm256_setzero_si256();
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m256i data_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
    __m256i mask = _mm256_cmpgt_epi32(data_vec, thresh_vec);
    __m256i ones = _mm256_and_si256(mask, _mm256_set1_epi32(1));
    count_vec = _mm256_add_epi32(count_vec, ones);
  }
  
  // Horizontal sum
  alignas(32) int32 counts[8];
  _mm256_store_si256(reinterpret_cast<__m256i*>(counts), count_vec);
  int32 total = 0;
  for (int i = 0; i < 8; ++i) {
    total += counts[i];
  }
  
  // Handle remaining elements
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] >= threshold) ++total;
  }
  
  return total;
#elif defined(HAS_SSE2)
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m128i thresh_vec = _mm_set1_epi32(threshold - 1);
  __m128i count_vec = _mm_setzero_si128();
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m128i data_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
    __m128i mask = _mm_cmpgt_epi32(data_vec, thresh_vec);
    __m128i ones = _mm_and_si128(mask, _mm_set1_epi32(1));
    count_vec = _mm_add_epi32(count_vec, ones);
  }
  
  // Horizontal sum
  alignas(16) int32 counts[4];
  _mm_store_si128(reinterpret_cast<__m128i*>(counts), count_vec);
  int32 total = counts[0] + counts[1] + counts[2] + counts[3];
  
  // Handle remaining elements
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] >= threshold) ++total;
  }
  
  return total;
#else
  // Fallback implementation
  int32 total = 0;
  for (int32 i = 0; i < size; ++i) {
    if (data[i] >= threshold) ++total;
  }
  return total;
#endif
}

// Vectorized increment of multiple indices
inline void increment_multiple_indices(int32* data, const std::vector<int32>& indices) {
#ifdef HAS_AVX2
  // For small numbers of indices, just do them individually
  // For larger numbers, we could use gather/scatter operations
  for (int32 idx : indices) {
    ++data[idx];
  }
#else
  for (int32 idx : indices) {
    ++data[idx];
  }
#endif
}

// SIMD-optimized find first element above threshold
inline int32 find_first_above_threshold(const int32* data, int32 size, int32 threshold) {
#ifdef HAS_AVX2
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m256i thresh_vec = _mm256_set1_epi32(threshold - 1);
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m256i data_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
    __m256i mask = _mm256_cmpgt_epi32(data_vec, thresh_vec);
    
    int mask_bits = _mm256_movemask_ps(_mm256_castsi256_ps(mask));
    if (mask_bits) {
      // Find first set bit
      int first_bit = __builtin_ctz(mask_bits);
      return i + first_bit;
    }
  }
  
  // Check remaining elements
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] >= threshold) return i;
  }
  
#elif defined(HAS_SSE2)
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m128i thresh_vec = _mm_set1_epi32(threshold - 1);
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m128i data_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
    __m128i mask = _mm_cmpgt_epi32(data_vec, thresh_vec);
    
    int mask_bits = _mm_movemask_ps(_mm_castsi128_ps(mask));
    if (mask_bits) {
      int first_bit = __builtin_ctz(mask_bits);
      return i + first_bit;
    }
  }
  
  // Check remaining elements
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] >= threshold) return i;
  }
  
#else
  for (int32 i = 0; i < size; ++i) {
    if (data[i] >= threshold) return i;
  }
#endif
  
  return -1; // Not found
}

// Batch comparison with multiple thresholds
inline void count_above_multiple_thresholds(const int32* data, int32 size, 
                                            const std::vector<int32>& thresholds,
                                            std::vector<int32>& results) {
  results.resize(thresholds.size());
  
  for (size_t t = 0; t < thresholds.size(); ++t) {
    results[t] = count_above_threshold(data, size, thresholds[t]);
  }
}

// SIMD horizontal maximum (find max element in array)
inline int32 find_maximum(const int32* data, int32 size) {
#ifdef HAS_AVX2
  if (size == 0) return 0;
  
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m256i max_vec = _mm256_set1_epi32(INT32_MIN);
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m256i data_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
    max_vec = _mm256_max_epi32(max_vec, data_vec);
  }
  
  // Horizontal maximum within the vector
  alignas(32) int32 maxes[8];
  _mm256_store_si256(reinterpret_cast<__m256i*>(maxes), max_vec);
  int32 result = maxes[0];
  for (int i = 1; i < 8; ++i) {
    if (maxes[i] > result) result = maxes[i];
  }
  
  // Check remaining elements
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] > result) result = data[i];
  }
  
  return result;
  
#elif defined(HAS_SSE2)
  if (size == 0) return 0;
  
  const int32 simd_size = (size / SIMD_WIDTH_32) * SIMD_WIDTH_32;
  __m128i max_vec = _mm_set1_epi32(INT32_MIN);
  
  for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
    __m128i data_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
    // SSE2 doesn't have max_epi32, so we need to emulate it
    __m128i mask = _mm_cmpgt_epi32(data_vec, max_vec);
    max_vec = _mm_or_si128(_mm_and_si128(mask, data_vec), _mm_andnot_si128(mask, max_vec));
  }
  
  alignas(16) int32 maxes[4];
  _mm_store_si128(reinterpret_cast<__m128i*>(maxes), max_vec);
  int32 result = maxes[0];
  for (int i = 1; i < 4; ++i) {
    if (maxes[i] > result) result = maxes[i];
  }
  
  for (int32 i = simd_size; i < size; ++i) {
    if (data[i] > result) result = data[i];
  }
  
  return result;
  
#else
  if (size == 0) return 0;
  int32 result = data[0];
  for (int32 i = 1; i < size; ++i) {
    if (data[i] > result) result = data[i];
  }
  return result;
#endif
}

} // namespace simd_utils

// Extended ClusterTable class with SIMD operations
class ClusterTableSIMD : public ClusterTable {
private:
  std::vector<int32> split_count_simd;  // SIMD-friendly storage
  // Align split_count for better SIMD performance
  alignas(32) std::vector<int32> split_count_simd_aligned;
  
  // Cache frequently accessed values
  mutable int32 cached_n_tip = -1;
  mutable int32 cached_threshold = -1;
  mutable int32 cached_above_threshold_count = -1;
  
public:
  // Explicit constructor that properly initializes the base class
  ClusterTableSIMD(const Rcpp::List& phylo) : ClusterTable(phylo) {
    // Base class should be fully constructed at this point
  }
  
  // Copy constructor
  ClusterTableSIMD(const ClusterTableSIMD& other) 
    : ClusterTable(static_cast<const ClusterTable&>(other)),
      split_count_simd(other.split_count_simd) {
  }
  
  // Assignment operator
  ClusterTableSIMD& operator=(const ClusterTableSIMD& other) {
    if (this != &other) {
      ClusterTable::operator=(other);
      split_count_simd = other.split_count_simd;
    }
    return *this;
  }
  
  // SIMD-specific methods
  inline void init_split_count_simd(int32 n_elements) {
    split_count_simd_aligned.clear();
    // Round up to SIMD width for better vectorization
    int32 aligned_size = ((n_elements + SIMD_WIDTH_32 - 1) / SIMD_WIDTH_32) *
      SIMD_WIDTH_32;
    split_count_simd_aligned.resize(aligned_size, 1);
    cached_n_tip = n_elements;
    cached_threshold = -1; // Invalidate cache
  }
  
  inline void increment_split_count_simd(int32 index) {
    if (index >= 0 && index < cached_n_tip) {
      ++split_count_simd_aligned[index];
      cached_threshold = -1; // Invalidate cache
    }
  }
  
  inline int32 get_split_count_simd(int32 index) const {
    if (index >= 0 && index < cached_n_tip) {
      return split_count_simd_aligned[index];
    }
    return 0;
  }
  
  inline int32 count_splits_above_threshold_simd(int32 threshold) const {
    if (split_count_simd.empty()) return 0;
    return simd_utils::count_above_threshold(
      split_count_simd.data(), 
      static_cast<int32>(split_count_simd.size()), 
      threshold
    );
  }
  // Add these methods to your ClusterTableSIMD class:
  
  // SIMD-optimized method to find all splits above threshold
  inline void find_splits_above_threshold_simd(int32 threshold, std::vector<int32>& result) const {
    result.clear();
    if (split_count_simd.empty()) return;
    
#ifdef HAS_AVX2
    const int32 simd_size = (split_count_simd.size() / SIMD_WIDTH_32) * SIMD_WIDTH_32;
    __m256i thresh_vec = _mm256_set1_epi32(threshold - 1);  // For >= comparison
    
    for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
      __m256i data_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&split_count_simd[i]));
      __m256i mask = _mm256_cmpgt_epi32(data_vec, thresh_vec);
      
      // Extract mask and check each element
      int mask_bits = _mm256_movemask_ps(_mm256_castsi256_ps(mask));
      for (int j = 0; j < SIMD_WIDTH_32; ++j) {
        if (mask_bits & (1 << j)) {
          result.push_back(i + j);
        }
      }
    }
    
    // Handle remaining elements
    for (int32 i = simd_size; i < static_cast<int32>(split_count_simd.size()); ++i) {
      if (split_count_simd[i] >= threshold) {
        result.push_back(i);
      }
    }
    
#elif defined(HAS_SSE2)
    const int32 simd_size = (split_count_simd.size() / SIMD_WIDTH_32) * SIMD_WIDTH_32;
    __m128i thresh_vec = _mm_set1_epi32(threshold - 1);
    
    for (int32 i = 0; i < simd_size; i += SIMD_WIDTH_32) {
      __m128i data_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&split_count_simd[i]));
      __m128i mask = _mm_cmpgt_epi32(data_vec, thresh_vec);
      
      int mask_bits = _mm_movemask_ps(_mm_castsi128_ps(mask));
      for (int j = 0; j < SIMD_WIDTH_32; ++j) {
        if (mask_bits & (1 << j)) {
          result.push_back(i + j);
        }
      }
    }
    
    // Handle remaining elements
    for (int32 i = simd_size; i < static_cast<int32>(split_count_simd.size()); ++i) {
      if (split_count_simd[i] >= threshold) {
        result.push_back(i);
      }
    }
    
#else
    // Fallback implementation
    for (int32 i = 0; i < static_cast<int32>(split_count_simd.size()); ++i) {
      if (split_count_simd[i] >= threshold) {
        result.push_back(i);
      }
    }
#endif
  }
  
  // Batch increment operations (for potential future optimization)
  inline void increment_split_counts_batch_simd(const std::vector<int32>& indices) {
    for (int32 idx : indices) {
      if (idx >= 0 && idx < static_cast<int32>(split_count_simd.size())) {
        ++split_count_simd[idx];
      }
    }
  }
  
  // SIMD-optimized bulk threshold check (returns count only)
  inline int32 count_splits_above_threshold_simd_bulk(int32 threshold) const {
    if (split_count_simd.empty()) return 0;
    return simd_utils::count_above_threshold(
      split_count_simd.data(), 
      static_cast<int32>(split_count_simd.size()), 
      threshold
    );
  }
  
  // Helper to get raw data pointer for SIMD operations
  inline const int32* get_split_count_data() const {
    return split_count_simd.empty() ? nullptr : split_count_simd.data();
  }
  
  inline int32 get_split_count_size() const {
    return static_cast<int32>(split_count_simd.size());
  }
  
  // Debug method to check if base class is properly initialized
  inline bool is_properly_initialized() const {
    try {
      // Cast away const to call non-const methods (safe for testing)
      ClusterTableSIMD* non_const_this = const_cast<ClusterTableSIMD*>(this);
      int32 n = non_const_this->N();  // Test base class method
      int32 m = non_const_this->M();  // Test another base class method
      return (n > 0 && m >= 0);
    } catch (...) {
      return false;
    }
  }
  
  // Cached threshold counting
  inline int32 count_splits_above_threshold_cached(int32 threshold) const {
    if (threshold != cached_threshold) {
      cached_threshold = threshold;
      cached_above_threshold_count = simd_utils::count_above_threshold(
        split_count_simd_aligned.data(), 
        cached_n_tip, 
        threshold
      );
    }
    return cached_above_threshold_count;
  }
  
  // Prefetch data for better cache performance
  inline void prefetch_split_data() const {
    if (!split_count_simd_aligned.empty()) {
#ifdef __builtin_prefetch
      const char* data_ptr = reinterpret_cast<const char*>(split_count_simd_aligned.data());
      int32 cache_lines = (cached_n_tip * sizeof(int32) + 63) / 64; // 64-byte cache lines
      for (int32 i = 0; i < cache_lines; ++i) {
        __builtin_prefetch(data_ptr + i * 64, 0, 3); // Read, high temporal locality
      }
#endif
    }
  }
};

} // namespace TreeTools

#endif // _TREETOOLS_CLUSTERTABLE_SIMD_H
