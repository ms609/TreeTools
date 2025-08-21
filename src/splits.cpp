#include <Rcpp/Lightest>
#include <memory> // for make_unique
#include <stdexcept> /* for errors */
#include "../inst/include/TreeTools/assert.h" /* for ASSERT */
#include "../inst/include/TreeTools.h"

using namespace Rcpp;

const uintx BIN_SIZE = 8;
const uintx bin_mask[BIN_SIZE + 1] = {255, 1, 3, 7, 15, 31, 63, 127, 255};

#define NOT_TRIVIAL TreeTools::UINTX_MAX

namespace TreeTools {
  template<typename T>
  constexpr void set_bit(T& target, int bit_pos) noexcept {
    target |= T(1) << bit_pos;
  }

  [[nodiscard]] constexpr uintx power_of_two(int bit_pos) noexcept {
    return uintx(1) << bit_pos;
  }
}
using TreeTools::power_of_two;

inline void check_16_bit(double x) {
  if (x > double(std::numeric_limits<int16>::max())) {
    Rcpp::stop("Cannot represent object this large in 16-bit register.");
  }
}

// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
Rcpp::RawMatrix cpp_edge_to_splits(const Rcpp::IntegerMatrix& edge,
                                   const Rcpp::IntegerVector& order,
                                   const Rcpp::IntegerVector& nTip) {
  
  // Check input is valid
  if (edge.cols() != 2) {
    Rcpp::stop("Edge matrix must contain two columns");
  }
  
  const uintx n_edge = edge.rows();
  if (n_edge + 1 >= NOT_TRIVIAL) {
    Rcpp::stop("Too many edges in tree for edge_to_splits: "                    // # nocov
                 "Contact maintainer for advice");                   // # nocov
  }
  
  if (nTip[0] < 1) {
    if (nTip[0] == 0) {
      return RawMatrix(0, 0);
    } else {
      Rcpp::stop("Tree must contain non-negative number of tips.");
    }
  }
  
  if (n_edge != static_cast<uintx>(order.length())) {
    Rcpp::stop("Length of `order` must equal number of edges");
  }
  
  const uintx n_node = n_edge + 1;
  const uintx n_tip = nTip[0];
  const uintx n_bin = ((n_tip - 1) / BIN_SIZE) + 1;
  
  if (n_edge == n_tip || n_tip < 4) {
    return RawMatrix(0, n_bin);
  }
  if (n_edge < 3) {
    /* Cannot calculate trivial_two below. */
    Rcpp::stop("Not enough edges in tree for edge_to_splits.");
  }
  
  std::vector<uintx> splits(n_node * n_bin, 0);
  
  auto split = [&](uintx i, uintx j) -> uintx& {
    return splits[i * n_bin + j];
  };
  
  // Tip initialization
  for (uintx i = 0; i < n_tip; ++i) {
    split(i, i / BIN_SIZE) = power_of_two(i % BIN_SIZE);
  }
  
  const int order_root = order[n_edge - 1];
  const uintx root_node = edge(order_root, 0);
  uintx root_child = edge(order_root, 1);
  int32 root_children = 1;
  
  for (uintx i = 0; i != n_edge - 1; ++i) { // Omit last edge
    const int order_i = order[i];
    const uintx parent = edge(order_i, 0);
    const uintx child  = edge(order_i, 1);
    
    if (parent == root_node) {
      ++root_children;
      if (child > n_tip) {
        root_child = child;
      }
    }
    
    uintx* parent_split = &splits[(parent - 1) * n_bin];
    const uintx* child_split = &splits[(child - 1) * n_bin];
    
    for (uintx j = 0; j < n_bin; ++j) {
      parent_split[j] |= child_split[j];
    }
  }
  
  const uintx trivial_origin = root_node - 1;
  const uintx trivial_two = (root_children == 2 ? root_child - 1 : NOT_TRIVIAL);
  const uintx n_return = n_edge - n_tip - (trivial_two != NOT_TRIVIAL ? 1 : 0);
  
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);
  
  std::vector<uintx> valid_rows;
  valid_rows.reserve(n_return);
  
  for (uintx i = n_tip; i < n_node; ++i) {
    if (i != trivial_origin && i != trivial_two) {
      valid_rows.push_back(i);
      names[valid_rows.size() - 1] = static_cast<int>(i + 1);
    }
  }
  
  Rbyte* __restrict__ ret_data = RAW(ret);
  const uintx* __restrict__ splits_data = splits.data();
  
  for (uintx j = 0; j < n_bin; ++j) {
    Rbyte* __restrict__ dest_col = ret_data + j * n_return;
    
    for (uintx r = 0; r < n_return; ++r) {
      dest_col[r] = static_cast<Rbyte>(splits_data[valid_rows[r] * n_bin + j]);
    }
  }
  
  rownames(ret) = names;
  return ret;
}

// [[Rcpp::export]]
LogicalVector duplicated_splits(const RawMatrix splits,
                                const LogicalVector fromLast) {
  if (!splits.hasAttribute("nTip")) {
    Rcpp::stop("`splits` lacks an `nTip` attribute.");
  }
  const IntegerVector nTip = as<IntegerVector>(splits.attr("nTip"));
  const intx
    n_split = splits.rows(),
    n_tip = nTip[0],
    n_bin = ((n_tip - 1) / BIN_SIZE) + 1,
    n_spare = n_tip % BIN_SIZE,
    check_bins = n_bin - (n_spare == 1)
  ;
  if (n_split == 0) {
    return Rcpp::LogicalVector(0);
  }
  // assert(n_split > 0);
  if (n_bin != splits.cols()) {
    Rcpp::stop("`splits` tip number is mis-specified.");
  }
  
  RawMatrix compare(n_split, check_bins);
  if (n_spare == 0) {
    for (intx i = n_split; i--; ) {
      if (splits(i, 0) & bin_mask[1]) {
        compare(i, _) = splits(i, _);
      } else {
        for (intx j = check_bins; j--; ) {
          compare(i, j) = ~splits(i, j);
        }
      }
    }
  } else if (n_spare == 1) {
    for (intx i = n_split; i--; ) {
      if (splits(i, n_bin - 1)) {
        for (intx j = check_bins; j--; ) {
          compare(i, j) = splits(i, j);
        }
      } else {
        for (intx j = check_bins; j--; ) {
          compare(i, j) = ~splits(i, j);
        }
      }
    }
  } else {
    for (intx i = n_split; i--; ) {
      if (splits(i, 0) % 2) {
        compare(i, check_bins - 1) = splits(i, check_bins - 1) ^
          decltype(splits(0, 0))(bin_mask[n_spare]);
        for (intx j = check_bins - 1; j--; ) {
          compare(i, j) = ~splits(i, j);
        }
      } else {
        compare(i, _) = splits(i, _);
      }
    }
  }
  
  LogicalVector ret(n_split);
  if (fromLast[0]) {
    for (intx it = n_split - 1; it--; ) {
      const intx i = it + 1; // nothing to duplicate split(0, _)
      if (ret[i]) {
        continue;
      }
      for (intx j = i; j--; ) {
        // Rcout << " check split " << i << " (" << uintx(compare(i, 0)) <<
        //   ") vs " << j << " (" << uintx(compare(j, 0)) << "): ";
        for(intx bin = 0; compare(i, bin) == compare(j, bin); ) {
          // Rcout << " [bin " << bin << "] ";
          ++bin;
          if (bin == check_bins) {
            // Rcout << "Duplicate!";
            ret[j] = true;
            break;
          }
        }
        // Rcout << "\n";
        
      }
    }
  } else {
    for (intx i = 0; i != n_split - 1; ++i) {
      if (ret[i]) {
        continue;
      }
      for (intx j = i + 1; j != n_split; ++j) {
        
        for(intx bin = 0; compare(i, bin) == compare(j, bin); ) {
          ++bin;
          if (bin == check_bins) {
            ret[j] = true;
            break;
          }
        }
        
      }
    }
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix mask_splits(RawMatrix x) {
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  check_16_bit(x.size());
  const int16
    n_tip = x.attr("nTip"),
    last_bin = int16(x.cols() - 1),
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  if (unset_tips == 0) {
    return x;
  } else {
    const uintx unset_mask = power_of_two(BIN_SIZE - unset_tips) - 1;
    for (int16 i = int16(x.rows()); i--; ) {
      x(i, last_bin) &= decltype(x(0, 0))(unset_mask);
    }
  }
  return x;
}

// [[Rcpp::export]]
RawMatrix not_splits(const RawMatrix x) {
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  check_16_bit(x.size());
  
  const int16
    n_tip = x.attr("nTip"),
    last_bin = int16(x.cols() - 1),
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  RawMatrix ret = clone(x);
  if (unset_tips == 0) {
    for (int16 i = int16(x.size()); i--; ) {
      ret[i] = ~ret[i];
    }
  } else {
    const uintx unset_mask = power_of_two(BIN_SIZE - unset_tips) - 1;
    for (int16 i = int16(x.rows()); i--; ) {
      ret(i, last_bin) = ~ret(i, last_bin) & decltype(ret(0, 0))(unset_mask);
    }
    for (int16 i = int16(x.rows()) * last_bin; i--; ) {
      ret[i] = ~ret[i];
    }
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix xor_splits(const RawMatrix x, const RawMatrix y) {
  check_16_bit(x.size());
  
  const int16 n_split = int16(x.rows());
  if (n_split != int16(y.rows())) {
    Rcpp::stop("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  const int16
    last_bin = int16(x.cols() - 1),
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  RawMatrix ret = clone(x);
  if (unset_tips == 0) {
    for (int16 i = int16(x.size()); i--; ) {
      ret[i] ^= y[i];
    }
  } else {
    const uintx unset_mask = power_of_two(BIN_SIZE - unset_tips) - 1;
    for (int16 i = int16(x.rows()); i--; ) {
      ret(i, last_bin) = (ret(i, last_bin) ^ y(i, last_bin)) & 
        decltype(ret(0, 0))(unset_mask);
    }
    for (int16 i = int16(x.rows()) * last_bin; i--; ) {
      ret[i] ^= y[i];
    }
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix and_splits(const RawMatrix x, const RawMatrix y) {
  check_16_bit(x.size());
  
  const int16 n_split = int16(x.rows());
  if (n_split != y.rows()) {
    Rcpp::stop("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  RawMatrix ret = clone(x);
  for (int16 i = int16(x.size()); i--; ) {
    ret[i] &= y[i];
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix or_splits(const RawMatrix x, const RawMatrix y) {
  check_16_bit(x.size());
  
  const int16 n_split = int16(x.rows());
  if (n_split != y.rows()) {
    Rcpp::stop("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  RawMatrix ret = clone(x);
  for (int16 i = int16(x.size()); i--; ) {
    ret[i] |= y[i];
  }
  return ret;
}

// [[Rcpp::export]]
Rcpp::List split_consistent(const RawMatrix needle,
                            const Rcpp::List haystacks,
                            const LogicalVector invert) {
  
  const bool CONTRADICTORY = invert[0] ? true : false;
  const bool CONSISTENT = invert[0] ? false : true;
  
  // Validate needle
  if (needle.rows() != 1) {
    Rcpp::stop("Needle must contain exactly one split (one row)");
  }
  if (!needle.hasAttribute("nTip")) {
    Rcpp::stop("Needle lacks nTip attribute");
  }
  
  const int16 n_tip = needle.attr("nTip");
  const int16 n_bin = ((n_tip - 1) / BIN_SIZE) + 1;
  
  if (needle.cols() != n_bin) {
    Rcpp::stop("Needle has incorrect number of columns for nTip"); // # nocov
  }
  
  const int n_haystacks = haystacks.size();
  Rcpp::List results(n_haystacks);
  
  // Process each haystack
  for (int h = 0; h < n_haystacks; h++) {
    if (TYPEOF(haystacks[h]) != RAWSXP) {
      Rcpp::stop("Haystack element %d is not a RawMatrix", h + 1);
    }
    
    RawMatrix haystack = haystacks[h];
    
    // Validate haystack
    if (!haystack.hasAttribute("nTip")) {
      Rcpp::stop("Haystack %d lacks nTip attribute", h + 1);
    }
    if (int16(haystack.attr("nTip")) != n_tip) {
      Rcpp::stop("Haystack %d has different nTip than needle", h + 1);
    }
    if (haystack.cols() != n_bin) {
      Rcpp::stop("Haystack %d has incorrect number of columns for nTip", h + 1); // # nocov
    }
    
    const int n_splits = haystack.rows();
    Rcpp::LogicalVector consistency(n_splits);
    
    // Check each split in the haystack against the needle
    for (int s = 0; s < n_splits; s++) {
      
      // Calculate the four intersections A∩C, A∩D, B∩C, B∩D
      // Where needle defines A|B and haystack[s] defines C|D
      // A = tips where needle bit = 1, B = tips where needle bit = 0
      // C = tips where haystack bit = 1, D = tips where haystack bit = 0
      
      bool ac_nonempty = false; // A ∩ C (both bits = 1)
      bool ad_nonempty = false; // A ∩ D (needle = 1, haystack = 0)
      bool bc_nonempty = false; // B ∩ C (needle = 0, haystack = 1)
      bool bd_nonempty = false; // B ∩ D (both bits = 0)
      
      // Check each tip by examining the appropriate bit
      for (int tip = 0; tip < n_tip; tip++) {
        const int bin_index = tip / 8;
        const uint8_t bit_mask = 1 << (tip % 8);
        
        const bool needle_bit = (needle(0, bin_index) & bit_mask) != 0;
        const bool haystack_bit = (haystack(s, bin_index) & bit_mask) != 0;
        
        // Update intersection flags based on bit combinations
        if (needle_bit && haystack_bit) {
          ac_nonempty = true;        // A ∩ C
        } else if (needle_bit && !haystack_bit) {
          ad_nonempty = true;        // A ∩ D
        } else if (!needle_bit && haystack_bit) {
          bc_nonempty = true;        // B ∩ C
        } else { // !needle_bit && !haystack_bit
          bd_nonempty = true;        // B ∩ D
        }
        
        // Early termination: if all four intersections are non-empty,
        // the splits are contradictory
        if (ac_nonempty && ad_nonempty && bc_nonempty && bd_nonempty) {
          break;
        }
      }
      
      // Two splits are contradictory if all four intersections are non-empty
      // Otherwise they are consistent
      if (ac_nonempty && ad_nonempty && bc_nonempty && bd_nonempty) {
        consistency[s] = CONTRADICTORY;
      } else {
        consistency[s] = CONSISTENT;
      }
    }
    
    results[h] = consistency;
  }
  
  return results;
}

// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix thin_splits(const RawMatrix splits, const LogicalVector drop) {
  if (static_cast<uintx>(drop.length()) > std::numeric_limits<uintx>::max()) {
    Rcpp::stop("Splits this large are not (yet) supported.");
  }
  
  // Initialize
  const uintx n_split = uintx(splits.rows()),
              n_tip = uintx(drop.length()),
              n_bin = ((n_tip - 1) / BIN_SIZE) + 1;
  ASSERT(int(n_bin) == splits.cols());

  
  RawMatrix ret(n_split, n_bin);
  auto in_split = std::make_unique<uintx[]>(n_split);
  uintx kept_tip = 0;
  // Populate splits
  for (uintx tip = 0; tip != n_tip; ++tip) {
    if (drop[tip]) {
      continue;
    }
    for (uintx split = n_split; split--; ) {
      if (splits(split, uintx(tip / BIN_SIZE))
            & power_of_two(tip % BIN_SIZE)) {
        TreeTools::set_bit(ret(split, uintx(kept_tip / BIN_SIZE)),
                           kept_tip % BIN_SIZE);
        ++in_split[split];
      }
    }
    ++kept_tip;
  }
  
  uintx kept_splits = n_split;
  for (uintx split = n_split; split--; ) {
    const uintx in_this = in_split[split];
    if (in_this < 2 || in_this > (kept_tip - 2)) {
      in_split[split] = 0;
      --kept_splits;
    }
  }

  if (!kept_tip) {
    return RawMatrix(kept_splits, 0);
  }
  
  const intx last_bin = (kept_tip - 1) / BIN_SIZE;
  if (!kept_splits) {
    return RawMatrix(0, last_bin + 1);
  }
  
  CharacterVector names = rownames(splits);
  if (last_bin) {
    ret = ret(_, Range(0, last_bin));
  } else {
    ret = RawMatrix(n_split, 1, ret(_, 0).begin());
  }
  
  if (kept_splits == n_split) {
    rownames(ret) = names;
    return ret;
  }
  
  if (kept_splits == 1) {
    for (uintx i = n_split; i--; ) {
      if (in_split[i]) {
        ret = RawMatrix(1, last_bin + 1, ret(i, _).begin());
        rownames(ret) = CharacterVector::create(names[i]);
        return ret;
      }
    }
  }
  
  CharacterVector new_names(kept_splits);
  uintx out = 0;
  for (uintx i = 0; out != kept_splits; ++i) {
    if (in_split[i]) {
      new_names[out] = names[i];
      ret(out, _) = ret(i, _);
      ++out;
    }
  }
  ret = ret(Range(0, out - 1), _);
  rownames(ret) = new_names;
  return ret;
  
}

// [[Rcpp::export]]
RawMatrix pack_splits_logical(LogicalMatrix x) {
  const int nrow = x.nrow();
  const int ncol = x.ncol();
  
  // Handle zero columns up-front (produces nrow x 0 matrix)
  if (ncol == 0) {
    return RawMatrix(nrow, 0);
  }
  
  const int nbin = (ncol + 7) / 8;  // ceil(ncol / 8)
  RawMatrix out(nrow, nbin);
  
  const int* xp = LOGICAL(x);       // column-major storage
  
  for (int r = 0; r < nrow; ++r) {
    for (int b = 0; b < nbin; ++b) {
      const int start_col = b * 8;
      int limit = ncol - start_col;
      if (limit > 8) limit = 8;
      
      // Pointer to (r, start_col) in column-major order
      const int* p = xp + r + start_col * nrow;
      
      unsigned char byte = 0;
      // Walk across columns in this bin by stepping one column = +nrow
      for (int k = 0; k < limit; ++k) {
        if (*p == TRUE) byte |= (1u << k);
        p += nrow;
      }
      out(r, b) = byte;
    }
  }
  return out;
}

// [[Rcpp::export]]
RawMatrix pack_splits_logical_vec(LogicalVector x) {
  const int ncol = x.size();
  
  // 1 x 0 for a length-0 vector
  if (ncol == 0) {
    return RawMatrix(1, 0);
  }
  
  const int nbin = (ncol + 7) / 8;
  RawMatrix out(1, nbin);
  
  const int* xp = LOGICAL(x);  // contiguous
  
  for (int b = 0; b < nbin; ++b) {
    const int start = b * 8;
    int limit = ncol - start;
    if (limit > 8) limit = 8;
    
    const int* p = xp + start;
    
    unsigned char byte = 0;
    for (int k = 0; k < limit; ++k) {
      if (p[k] == TRUE) byte |= (1u << k);
    }
    out(0, b) = byte;
  }
  return out;
}

// Fast count of unique splits in a phylogenetic tree
// [[Rcpp::export]]
int cpp_count_splits(const Rcpp::IntegerMatrix& edge, const int nTip) {
  if (edge.ncol() != 2) {
    Rcpp::stop("Edge matrix must contain two columns");
  }
  
  const int n_edge = edge.nrow();
  if (n_edge <= nTip || nTip < 4) {
    return 0;
  }
  
  // Find max node number to size our vectors
  int max_node = nTip;  // Start with nTip as minimum
  for (int i = 0; i < n_edge; ++i) {
    max_node = std::max(max_node, std::max(edge(i, 0), edge(i, 1)));
  }
  
  // Use vectors instead of maps - much faster O(1) access
  std::vector<int> parent_count(max_node + 1, 0);
  std::vector<int> child_count(max_node + 1, 0);
  std::vector<bool> is_internal(max_node + 1, false);
  
  // Single pass to count occurrences and mark internal nodes
  for (int i = 0; i < n_edge; ++i) {
    const int parent = edge(i, 0);
    const int child = edge(i, 1);
    
    parent_count[parent]++;
    child_count[child]++;
    
    // Mark internal nodes (> nTip)
    if (parent > nTip) {
      is_internal[parent] = true;
    }
    if (child > nTip) {
      is_internal[child] = true;
    }
  }
  
  // Count internal nodes and singletons in single pass
  int n_internal = 0;
  int n_singles = 0;
  int root_node = 0;
  
  for (int node = nTip + 1; node <= max_node; ++node) {
    if (is_internal[node]) {
      n_internal++;
      
      // Check if singleton (degree 2: appears once as parent AND once as child)
      if (parent_count[node] == 1 && child_count[node] == 1) {
        n_singles++;
      }
      
      // Check if root (appears as parent but not as child)
      if (child_count[node] == 0) {
        root_node = node;
      }
    }
  }
  
  if (n_internal == 0) {
    return 0;
  }
  
  // TreeIsRooted: root has < 3 children
  const bool is_rooted = parent_count[root_node] < 3;
  
  return (n_internal - n_singles) - 1 - (is_rooted ? 1 : 0);
}
