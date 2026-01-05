#ifndef _TREETOOLS_SPLITLIST_H
#define _TREETOOLS_SPLITLIST_H

#include <Rcpp/Lightest>
#include <stdexcept> /* for errors */
#include <vector>    /* for heap allocation */
#include <algorithm> /* for std::fill */

#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16, int32 */

using splitbit = uint_fast64_t;

#define R_BIN_SIZE int16(8)
#define SL_BIN_SIZE int16(64)
#define SL_MAX_BINS int16(32)

/* * Stack allocation limits (Legacy support for speed)
 * Trees smaller than this will use stack arrays.
 * Trees larger will trigger heap allocation.
 */
#define SL_MAX_TIPS (SL_BIN_SIZE * SL_MAX_BINS) // 2048
#define SL_MAX_SPLITS (SL_MAX_TIPS - 3) 

#define INLASTBIN(n, size) int16((size) - int16((size) - int16((n) % (size))) % (size))
#define INSUBBIN(bin, offset)                                  \
  splitbit(x(split, ((bin) * input_bins_per_bin) + (offset)))
#define INBIN(r_bin, bin) ((INSUBBIN((bin), (r_bin))) << (R_BIN_SIZE * (r_bin)))

namespace TreeTools {

  constexpr int input_bins_per_bin = SL_BIN_SIZE / R_BIN_SIZE;

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
    int32 n_splits;
    int16 n_bins;
    int16* in_split;
    splitbit** state;
  
  private:
    /* * STACK STORAGE (Fast path for small trees) */
    int16 stack_in_split[SL_MAX_SPLITS];
    splitbit stack_state[SL_MAX_SPLITS][SL_MAX_BINS];
    splitbit* stack_rows[SL_MAX_SPLITS]; 
    
    /* * HEAP STORAGE (Large trees) */
    std::vector<int16> heap_in_split;
    std::vector<splitbit> heap_data;      
    std::vector<splitbit*> heap_rows;     
    
  public:
    SplitList(const Rcpp::RawMatrix &x) {
      
      const double n_rows = static_cast<double>(x.rows());
      
      /* Check limits */
      if (n_rows > static_cast<double>(std::numeric_limits<int32>::max())) {
        Rcpp::stop("Too many splits (exceeds int32 limit).");
      }
      
      n_splits = int32(x.rows());
      
      // Handle empty case safely
      if (n_splits == 0) {
        n_bins = 0;
        in_split = nullptr; // or stack_in_split
        state = nullptr;    // or stack_rows
        return;
      }
      
      const int16 n_input_bins = int16(x.cols());
      
      if (n_input_bins == 0) {
        n_bins = 0;
      } else {
        n_bins = (n_input_bins + input_bins_per_bin - 1) / input_bins_per_bin;
      }
      
      bool use_heap = (n_splits > SL_MAX_SPLITS) || (n_bins > SL_MAX_BINS);
      
      if (use_heap) {
        heap_in_split.resize(n_splits, 0);
        in_split = heap_in_split.data();
        
        size_t total_elements = static_cast<size_t>(n_splits) * static_cast<size_t>(n_bins);
        heap_data.resize(total_elements); 
        
        heap_rows.resize(n_splits);
        for (int32 i = 0; i < n_splits; ++i) {
          heap_rows[i] = &heap_data[i * n_bins];
        }
        state = heap_rows.data();
        
      } else {
        in_split = stack_in_split;
        for (int16 i = 0; i < n_splits; ++i) {
          stack_rows[i] = stack_state[i];
          in_split[i] = 0; 
        }
        state = stack_rows;
      }
      
      // If no bins (no tips/data), we are done
      if (n_bins == 0) return;
      
      // 1. Process full bins (if any)
      // We iterate up to n_bins - 1. If n_bins is 1, this loop doesn't run.
      for (int16 bin = 0; bin < n_bins - 1; ++bin) {
        const int16 bin_offset = bin * input_bins_per_bin;
        
        for (int32 split = 0; split < n_splits; ++split) {
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
      const int16 bin_offset = last_bin * input_bins_per_bin;
      
      const int16 items_in_last = n_input_bins - bin_offset;
      
      for (int32 split = 0; split < n_splits; ++split) {
        splitbit combined = splitbit(x(split, bin_offset));
        
        for (int16 i = 1; i < items_in_last; ++i) {
          combined |= splitbit(x(split, bin_offset + i)) << (R_BIN_SIZE * i);
        }
        
        state[split][last_bin] = combined;
        in_split[split] += count_bits(combined);
      }
    }
    
    // Default destructor handles vector cleanup automatically.
    ~SplitList() = default;
    
    // Disable copy/move to prevent pointer invalidation issues 
    SplitList(const SplitList&) = delete;
    SplitList& operator=(const SplitList&) = delete;
  };
}

#endif
