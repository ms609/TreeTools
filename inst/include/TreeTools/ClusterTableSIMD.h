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

} // namespace simd_utils

// Extended ClusterTable class with SIMD operations
class ClusterTableSIMD : public ClusterTable {
private:
  std::vector<int32> split_count_simd;  // SIMD-friendly storage
  
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
    split_count_simd.clear();
    split_count_simd.resize(n_elements, 1);  // Initialize to 1 like original
  }
  
  inline void increment_split_count_simd(int32 index) {
    if (index >= 0 && index < static_cast<int32>(split_count_simd.size())) {
      ++split_count_simd[index];
    }
  }
  
  inline int32 get_split_count_simd(int32 index) const {
    if (index >= 0 && index < static_cast<int32>(split_count_simd.size())) {
      return split_count_simd[index];
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
};

} // namespace TreeTools

#endif // _TREETOOLS_CLUSTERTABLE_SIMD_H