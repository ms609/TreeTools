#ifndef _TREETOOLS_SPLITLIST_H
#define _TREETOOLS_SPLITLIST_H

#include <Rcpp/Lightest>

#include <stdexcept> /* for errors */

#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16 */

#define R_BIN_SIZE int16(8)
#define SL_BIN_SIZE int16(64)
#define SL_MAX_BINS int16(32)
/* 64*32 is about the largest size for which two SplitList objects reliably fit
 * on the stack (as required in TreeDist; supporting more leaves would mean
 * refactoring to run on the heap (and, trivially, converting int16 to int32
 * for split*bin implicit calculation in state[split][bin]?) */
#define SL_MAX_TIPS (SL_BIN_SIZE * SL_MAX_BINS)
#define SL_MAX_SPLITS (SL_MAX_TIPS - 3) /* no slower than a power of two */

#define INLASTBIN(n, size) int16((size) - int16((size) - int16((n) % (size))) % (size))
#define INSUBBIN(bin, offset)                                  \
  splitbit(x(split, ((bin) * input_bins_per_bin) + (offset)))
#define INBIN(r_bin, bin) ((INSUBBIN((bin), (r_bin))) << (R_BIN_SIZE * (r_bin)))

#define right16bits splitbit(65535U)

#define TREETOOLS_SPLITLIST_INIT __attribute__((constructor))  \
  void _treetools_initialize_bitcounts() {                     \
  for (int i = 65536; i--; ) {                                 \
    int16 n_bits = 0;                                          \
    for (int j = 16; j--; ) {                                  \
      if (i & (1 << j)) n_bits += 1;                           \
    }                                                          \
    TreeTools::bitcounts[i] = n_bits;                          \
  }                                                            \
}                                                              \

typedef uint_fast64_t splitbit;

namespace TreeTools {

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

  // Static here means that each translation unit (i.e. file, resulting in an .o)
  // will have its own copy of the variable (which it will initialize separately).
  static int16 bitcounts[65536];

  inline int16 count_bits (splitbit x) {
    /* For 32-bit splitbits: */
    /* return bitcounts[x & right16bits] + bitcounts[x >> 16]; */

    /* For 64-bit splitbits: */
    return int16(bitcounts[x & right16bits] + bitcounts[(x >> 16) & right16bits]
    + bitcounts[(x >> 32) & right16bits]
    + bitcounts[(x >> 48)]);
  }


  class SplitList {
  public:
    int16 n_splits, n_bins;
    int16 in_split[SL_MAX_SPLITS];
    splitbit state[SL_MAX_SPLITS][SL_MAX_BINS];
    SplitList(Rcpp::RawMatrix x) {
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
      
      const int16 n_input_bins = int16(x.cols()),
        input_bins_per_bin = SL_BIN_SIZE / R_BIN_SIZE;

      n_bins = int16(n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;

      if (n_bins > SL_MAX_BINS) {
        Rcpp::stop("This many leaves cannot be supported. "
                   "Please contact the TreeTools maintainer if "
                   "you need to use more!");
      }

      for (int16 split = 0; split != n_splits; split++) {
        int16 last_bin = n_bins - 1;
        const int16 raggedy_bins = INLASTBIN(n_input_bins, R_BIN_SIZE);
        /*Rcout << n_input_bins << " bins in; " << raggedy_bins << " raggedy bins\n";*/
        state[split][last_bin] = INSUBBIN(last_bin, 0);
        /*Rcout << " State[" << split << "][" << bin << "] = " << state[split][bin] << ".\n";*/
        for (int16 input_bin = 1; input_bin != raggedy_bins; input_bin++) {
          /*Rcout << "Adding " << (splitbit (x(split, (bin * input_bins_per_bin) + input_bin))) << " << "
                  << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" << bin
                  << "], was " << state[split][bin] << "\n";*/
          state[split][last_bin] += INBIN(input_bin, last_bin);
        }
        in_split[split] = count_bits(state[split][last_bin]);

        for (int16 bin = 0; bin != n_bins - 1; bin++) {
          /*Rcout << "Split " << split << ", bin << " << bin << ".\n";*/
          state[split][bin] = INSUBBIN(bin, 0);
          for (int16 input_bin = 1; input_bin != input_bins_per_bin; input_bin++) {
            /*Rcout << "Adding " << (splitbit (x(split, (bin * input_bins_per_bin) + input_bin))) << " << "
                    << (R_BIN_SIZE * input_bin) << " to state [" << split << "]["
                    << bin << "], was " << state[split][bin] << "\n";*/
            state[split][bin] += INBIN(input_bin, bin);
          }
          in_split[split] += count_bits(state[split][bin]);
        }
      }
    }
  };
}

#endif
