#ifndef _TREETOOLS_SPLITLIST_H
#define _TREETOOLS_SPLITLIST_H

#include <Rcpp/Lightest>

#include <stdexcept> /* for errors */

#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16 */

using splitbit = uint_fast64_t;

#define R_BIN_SIZE int16(8)
#define SL_BIN_SIZE int16(64)
#define SL_MAX_BINS int16(32)
/* 64*32 is about the largest size for which two SplitList objects reliably fit
 * on the stack (as required in TreeDist; supporting more leaves would mean
 * refactoring to run on the heap (and, trivially, converting int16 to int32
 * for split*bin implicit calculation in state[split][bin]?) */
#define SL_MAX_TIPS (SL_BIN_SIZE * SL_MAX_BINS) // 32 * 64 = 2048
#define SL_MAX_SPLITS (SL_MAX_TIPS - 3) /* no slower than a power of two */

#define INLASTBIN(n, size) int16((size) - int16((size) - int16((n) % (size))) % (size))
#define INSUBBIN(bin, offset)                                  \
  splitbit(x(split, ((bin) * input_bins_per_bin) + (offset)))
#define INBIN(r_bin, bin) ((INSUBBIN((bin), (r_bin))) << (R_BIN_SIZE * (r_bin)))

// Retained for backward compatibility; not required since 1.15.0.9006
#define TREETOOLS_SPLITLIST_INIT __attribute__((constructor))  \
void _treetools_initialize_bitcounts() {                       \
  for (int i = 65536; i--; ) {                                 \
    int16 n_bits = 0;                                          \
    for (int j = 16; j--; ) {                                  \
      if (i & (1 << j)) n_bits += 1;                           \
    }                                                          \
    TreeTools::bitcounts[i] = n_bits;                          \
  }                                                            \
}                                                              \
  
namespace TreeTools {

  constexpr int input_bins_per_bin = SL_BIN_SIZE / R_BIN_SIZE;

// Retained for backward compatibility; not required since 1.15.0.9006
  static int16 bitcounts[65536];
  
// Retained for backward compatibility; not required since 1.15.0.9006
  const splitbit powers_of_two[SL_BIN_SIZE] = {
    0x1, 0x2, 0x4, 0x8,
    0x10, 0x20, 0x40, 0x80,
    0x100, 0x200, 0x400, 0x800,
    0x1000, 0x2000, 0x4000, 0x8000,
    0x10000, 0x20000, 0x40000, 0x80000,
    0x100000, 0x200000, 0x400000, 0x800000,
    0x1000000, 0x2000000, 0x4000000, 0x8000000,
    0x10000000, 0x20000000, 0x40000000, 0x80000000,
    0x100000000, 0x200000000, 0x400000000, 0x800000000,
    0x1000000000, 0x2000000000, 0x4000000000, 0x8000000000,
    0x10000000000, 0x20000000000, 0x40000000000, 0x80000000000,
    0x100000000000, 0x200000000000, 0x400000000000, 0x800000000000,
    0x1000000000000, 0x2000000000000, 0x4000000000000, 0x8000000000000,
    0x10000000000000, 0x20000000000000, 0x40000000000000, 0x80000000000000,
    0x100000000000000, 0x200000000000000, 0x400000000000000, 0x800000000000000,
    0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000
  };
  template<typename T>
  [[nodiscard]] constexpr T power_of_two(int bit_pos) noexcept {
    static_assert(std::is_unsigned_v<T>, "Use unsigned types for bit operations");
    assert(bit_pos >= 0 && bit_pos < int(sizeof(T) * 8));
    return T(1) << bit_pos;
  }
  
#if __cplusplus >= 202002L
#include <bit> // C++20 header for std::popcount
  
  inline int16 count_bits(splitbit x) {
    return static_cast<int16>(std::popcount(x));
  }
  
  // Option 2: Fallback for C++17 and older
#else
#if defined(__GNUC__) || defined(__clang__)
  // GCC and Clang support __builtin_popcountll for long long
  inline int16 count_bits(splitbit x) {
    return static_cast<int16>(__builtin_popcountll(x));
  }
#elif defined(_MSC_VER)
#include <intrin.h>
  inline int16 count_bits(splitbit x) {
    return static_cast<int16>(__popcnt64(x));
  }
#else
  // A slower, but safe and highly portable fallback for all other compilers
  // This is a last resort if no built-in is available.
  inline int16_t count_bits(splitbit x) {
    int16_t count = 0;
    while (x != 0) {
      x &= (x - 1);
      count++;
    }
    return count;
  }
#endif // Compiler check for builtins
  
#endif // C++20 check

  class SplitList {
  public:
    int16 n_splits, n_bins;
    int16 in_split[SL_MAX_SPLITS];
    splitbit state[SL_MAX_SPLITS][SL_MAX_BINS];
    SplitList(const Rcpp::RawMatrix &x) {
      if (double(x.rows()) > double(std::numeric_limits<int16>::max())) {
        Rcpp::stop("This many splits cannot be supported. "
                   "Please contact the TreeTools maintainer if "
                   "you need to use more!");
      }
      if (double(x.cols()) > double(std::numeric_limits<int16>::max())) {
        Rcpp::stop("This many leaves cannot be supported. "
                     "Please contact the TreeTools maintainer if "
                     "you need to use more!");
      }
      
      n_splits = int16(x.rows());
      ASSERT(n_splits >= 0);
      
      const int16 n_input_bins = int16(x.cols());

      n_bins = int16(n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;

      if (n_bins > SL_MAX_BINS) {
        Rcpp::stop("This many leaves cannot be supported. "
                   "Please contact the TreeTools maintainer if "
                   "you need to use more!");
      }
      
      for (int16 split = 0; split < n_splits; ++split) {
        in_split[split] = 0;
      }
      
      for (int16 bin = 0; bin < n_bins - 1; ++bin) {
        const int16 bin_offset = bin * input_bins_per_bin;
        
        for (int16 split = 0; split < n_splits; ++split) {
          splitbit combined = splitbit(x(split, bin_offset));
          
          for (int16 input_bin = 1; input_bin < input_bins_per_bin; ++input_bin) {
            combined |= splitbit(x(split, bin_offset + input_bin)) <<
              (R_BIN_SIZE * input_bin);
          }
          
          state[split][bin] = combined;
          in_split[split] += count_bits(combined);
        }
      }
      
      const int16 last_bin = n_bins - 1;
      const int16 raggedy_bins = INLASTBIN(n_input_bins, R_BIN_SIZE);
      
      for (int16 split = 0; split < n_splits; ++split) {
        state[split][last_bin] = INSUBBIN(last_bin, 0);
        
        for (int16 input_bin = 1; input_bin < raggedy_bins; ++input_bin) {
          state[split][last_bin] += INBIN(input_bin, last_bin);
        }
        
        in_split[split] += count_bits(state[split][last_bin]);
      }
    }
  };
}

#endif
