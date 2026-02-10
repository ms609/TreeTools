#include <Rcpp.h>
#include <unordered_map>
#include <string>

using namespace Rcpp;

// Helper: fill `key` with the canonicalized bytes for row `i`
// so that complements map to the same key. We follow the logic pattern
// in duplicated_splits(): choose an orientation based on a sentinel bit,
// and when inverting, fix the spare bits in the last used bin.
static inline void canonicalize_row_key(std::string &key,
                                        const RawMatrix &M,
                                        int i,
                                        int check_bins,
                                        int n_bin,
                                        int n_spare) {
  key.resize(check_bins);
  
  // mask with ones in the UNUSED bits of the last bin, when n_spare > 0.
  unsigned char unused_mask = 0u;
  if (n_spare > 0) {
    unsigned int used_mask = (1u << n_spare) - 1u;        // low n_spare bits set
    unused_mask = static_cast<unsigned char>(~used_mask); // high bits set
  }
  
  if (n_spare == 0) {
    // Decide orientation by LSB of first bin
    const bool keep = (static_cast<unsigned char>(M(i, 0)) & 0x01u) != 0u;
    if (keep) {
      for (int b = 0; b < check_bins; ++b) key[b] = static_cast<char>(M(i, b));
    } else {
      for (int b = 0; b < check_bins; ++b)
        key[b] = static_cast<char>(~static_cast<unsigned char>(M(i, b)));
    }
  } else if (n_spare == 1) {
    // Decide orientation by non-zero last bin
    const bool keep = static_cast<unsigned char>(M(i, n_bin - 1)) != 0u;
    if (keep) {
      for (int b = 0; b < check_bins; ++b) key[b] = static_cast<char>(M(i, b));
    } else {
      for (int b = 0; b < check_bins; ++b)
        key[b] = static_cast<char>(~static_cast<unsigned char>(M(i, b)));
    }
  } else {
    // Multiple spare bits:
    // If LSB of first bin is 1, invert all bins; after inversion,
    // fix the last bin's unused bits by XOR with `unused_mask`.
    const bool invert = (static_cast<unsigned char>(M(i, 0)) & 0x01u) != 0u;
    if (!invert) {
      for (int b = 0; b < check_bins; ++b) key[b] = static_cast<char>(M(i, b));
    } else {
      for (int b = 0; b < check_bins - 1; ++b)
        key[b] = static_cast<char>(~static_cast<unsigned char>(M(i, b)));
      unsigned char last = static_cast<unsigned char>(~static_cast<unsigned char>(M(i, check_bins - 1)));
      last ^= unused_mask; // zero-out the unused bits after inversion
      key[check_bins - 1] = static_cast<char>(last);
    }
  }
}

// [[Rcpp::export]]
IntegerVector first_matching_split_pair(const RawMatrix x, const RawMatrix table) {
  // Validate attributes and shapes
  if (!x.hasAttribute("nTip") || !table.hasAttribute("nTip"))
    stop("Both `x` and `table` must have an `nTip` attribute.");
  
  const int n_tip_x = as<IntegerVector>(x.attr("nTip"))[0];
  const int n_tip_t = as<IntegerVector>(table.attr("nTip"))[0];
  if (n_tip_x != n_tip_t)
    stop("`x` and `table` must have the same number of tips.");
  
  const int n_split_x = x.rows();
  const int n_split_t = table.rows();
  if (n_split_x == 0 || n_split_t == 0) return IntegerVector::create(0, 0);
  
  const int n_bin_x = x.cols();
  const int n_bin_t = table.cols();
  if (n_bin_x != n_bin_t)
    stop("`x` and `table` have incompatible bin counts.");
  const int n_bin = n_bin_x;
  
  // Compute bin layout
  const int BIN_SIZE = 8;
  const int expected_n_bin = ((n_tip_x - 1) / BIN_SIZE) + 1;
  if (expected_n_bin != n_bin)
    stop("`nTip` inconsistent with number of bins.");
  
  const int n_spare = n_tip_x % BIN_SIZE;
  const int check_bins = n_bin - ((n_spare == 1) ? 1 : 0);
  
  // Build hash from table
  std::unordered_map<std::string, int> H;
  H.reserve(static_cast<size_t>(n_split_t) * 2u);
  
  std::string key;
  key.reserve(static_cast<size_t>(check_bins));
  
  for (int j = 0; j < n_split_t; ++j) {
    canonicalize_row_key(key, table, j, check_bins, n_bin, n_spare);
    // keep first occurrence
    H.emplace(key, j + 1);
  }
  
  // Probe with x; return first hit (1-based indices)
  for (int i = 0; i < n_split_x; ++i) {
    canonicalize_row_key(key, x, i, check_bins, n_bin, n_spare);
    auto it = H.find(key);
    if (it != H.end()) {
      return IntegerVector::create(i + 1, it->second);
    }
  }
  
  return IntegerVector::create(0, 0);
}

// Convenience: only the index in `x` (0 if none)
// [[Rcpp::export]]
int first_matching_split_index(const RawMatrix x, const RawMatrix table) {
  IntegerVector ij = first_matching_split_pair(x, table);
  return ij[0];
}
