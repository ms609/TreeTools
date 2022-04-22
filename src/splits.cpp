#include <Rcpp/Lightest>
#include <memory> // for make_unique
#include <stdexcept> /* for errors */
#include "../inst/include/TreeTools/assert.h"
#include "../inst/include/TreeTools.h"
using namespace Rcpp;

const uintx powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                                 2048, 4096, 8192, 16384, 32768};
const uintx BIN_SIZE = 8;
const uintx bin_mask[BIN_SIZE + 1] = {255, 1, 3, 7, 15, 31, 63, 127, 255};

#define NOT_TRIVIAL TreeTools::UINTX_MAX

#define PO_PARENT(i) edge(order[i], 0)
#define PO_CHILD(i) edge(order[i], 1)

// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix cpp_edge_to_splits(const IntegerMatrix edge,
                             const IntegerVector order,
                             const IntegerVector nTip) {
  // Check input is valid
  if (edge.cols() != 2) {
    Rcpp::stop("Edge matrix must contain two columns");
  }
  if (1UL + edge.rows() > NOT_TRIVIAL - 1U) {
    Rcpp::stop("Too many edges in tree for edge_to_splits: "       // # nocov
                            "Contact maintainer for advice");                   // # nocov
  }
  if (nTip[0] < 1) {
    Rcpp::stop("Tree must contain tips.");
  }
  
  const uintx n_edge = edge.rows();
  if (n_edge != uintx(order.length())) {
    Rcpp::stop("Length of `order` must equal number of edges");
  }
  
  // Initialize
  const uintx n_node = n_edge + 1,
              root_node = PO_PARENT(n_edge - 1),
              n_tip = nTip[0],
              n_bin = ((n_tip - 1) / BIN_SIZE) + 1;

  if (n_edge == n_tip  /* No internal nodes resolved */
        || n_tip < 4) { /* Need four tips to split non-trivially */
    return RawMatrix (0, n_bin);
  }
  if (n_edge < 3) {
    /* Cannot calculate trivial_two below. */
    Rcpp::stop("Not enough edges in tree for edge_to_splits.");
  }

  uintx** splits = new uintx*[n_node];
  for (uintx i = n_node; i--; ) {
    splits[i] = new uintx[n_bin](); // () zero-initializes
  }

  // Populate splits
  for (uintx i = n_tip; i--; ) {
    splits[i][uintx(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }

  uintx root_child = PO_CHILD(n_edge - 1);
  int32 root_children = 1;
  for (uintx i = 0; i != n_edge - 1; ++i) { // Omit last edge
    const uintx parent = PO_PARENT(i);
    const uintx child = PO_CHILD(i);
    if (parent == root_node) {
      ++root_children;
      if (child > n_tip) {
        root_child = uintx(child);
      }
    }
    for (uintx j = 0; j != n_bin; ++j) {
      splits[parent - 1][j] |= splits[child - 1][j];
    }
  }

  for (uintx i = n_tip; i --; ) {
    delete[] splits[i];
  }

  // Only return non-trivial splits
  uintx n_trivial = 0;
  const uintx trivial_origin = root_node - 1,
              trivial_two = (root_children == 2 ? root_child - 1 : NOT_TRIVIAL),
              n_return = n_edge - n_tip - (trivial_two != NOT_TRIVIAL ? 1 : 0);
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);

  for (uintx i = n_tip; i != n_node; ++i) {
    if (i == trivial_origin || i == trivial_two) {
      ++n_trivial;
    } else {
      for (uintx j = 0; j != n_bin; j++) {
        ret(i - n_tip - n_trivial, j) = splits[i][j];
        names[i - n_tip - n_trivial] = (i + 1);
      }
    }
    delete[] splits[i];
  }

  delete[] splits;

  rownames(ret) = names;
  return(ret);
}

// [[Rcpp::export]]
LogicalVector duplicated_splits(const RawMatrix splits,
                                const LogicalVector fromLast) {
  if (!splits.hasAttribute("nTip")) {
    Rcpp::stop("`splits` lacks an `nTip` attribute.");
  }
  const IntegerVector nTip = splits.attr("nTip");
  const intx
    n_split = splits.rows(),
    n_tip = nTip[0],
    n_bin = ((n_tip - 1) / BIN_SIZE) + 1,
    n_spare = n_tip % BIN_SIZE,
    check_bins = n_bin - (n_spare == 1)
  ;
  if (n_bin != splits.cols()) {
    Rcpp::stop("`splits` tip number is mis-specified.");
  }
  if (n_split == 0) {
    return Rcpp::LogicalVector(0);
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
          bin_mask[n_spare];
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
  const int16
    n_tip = x.attr("nTip"),
    last_bin = x.cols() - 1,
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  if (unset_tips == 0) {
    return x;
  } else {
    const uintx unset_mask = powers_of_two[BIN_SIZE - unset_tips] - 1;
    for (int16 i = x.rows(); i--; ) {
      x(i, last_bin) &= unset_mask;
    }
  }
  return x;
}

// [[Rcpp::export]]
RawMatrix not_splits(const RawMatrix x) {
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  const int16
    n_tip = x.attr("nTip"),
    last_bin = x.cols() - 1,
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  RawMatrix ret = clone(x);
  if (unset_tips == 0) {
    for (int16 i = x.size(); i--; ) {
      ret[i] = ~ret[i];
    }
  } else {
    const uintx unset_mask = powers_of_two[BIN_SIZE - unset_tips] - 1;
    for (int16 i = x.rows(); i--; ) {
      ret(i, last_bin) = ~ret(i, last_bin) & unset_mask;
    }
    for (int16 i = x.rows() * last_bin; i--; ) {
      ret[i] = ~ret[i];
    }
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix xor_splits(const RawMatrix x, const RawMatrix y) {
  const int16 n_split = x.rows();
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
  
  const int16
    last_bin = x.cols() - 1,
    unset_tips = (n_tip % BIN_SIZE) ? BIN_SIZE - n_tip % BIN_SIZE : 0
  ;
  
  RawMatrix ret = clone(x);
  if (unset_tips == 0) {
    for (int16 i = x.size(); i--; ) {
      ret[i] ^= y[i];
    }
  } else {
    const uintx unset_mask = powers_of_two[BIN_SIZE - unset_tips] - 1;
    for (int16 i = x.rows(); i--; ) {
      ret(i, last_bin) = (ret(i, last_bin) ^ y(i, last_bin)) & unset_mask;
    }
    for (int16 i = x.rows() * last_bin; i--; ) {
      ret[i] ^= y[i];
    }
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix and_splits(const RawMatrix x, const RawMatrix y) {
  const int16 n_split = x.rows();
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
  for (int16 i = x.size(); i--; ) {
    ret[i] &= y[i];
  }
  return ret;
}

// [[Rcpp::export]]
RawMatrix or_splits(const RawMatrix x, const RawMatrix y) {
  const int16 n_split = x.rows();
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
  for (int16 i = x.size(); i--; ) {
    ret[i] |= y[i];
  }
  return ret;
}

// Edges must be listed in 'strict' postorder, i.e. two-by-two
// [[Rcpp::export]]
RawMatrix thin_splits(const RawMatrix splits, const LogicalVector drop) {
    // Initialize
  const uintx n_split = splits.rows(),
              n_tip = drop.length(),
              n_bin = ((n_tip - 1) / BIN_SIZE) + 1;
  ASSERT(n_bin == splits.cols());

  
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
            & powers_of_two[tip % BIN_SIZE]) {
        ret(split, uintx(kept_tip / BIN_SIZE)) |= 
          powers_of_two[kept_tip % BIN_SIZE];
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
