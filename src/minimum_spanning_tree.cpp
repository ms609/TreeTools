#include <Rcpp.h>
#include <memory> /* for unique_ptr */
#include <cmath> /* for ceil, sqrt() */
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;
using namespace std;

intx island_housing(intx x, unique_ptr<intx[]> &island) {
  if (x == island[x]) return x;
  intx parent_house = island_housing(island[x], island);
  if (parent_house == island[x]) return island[x];
  island[x] = parent_house;
  return parent_house;
}

// @param order Index of largest, second largest, ... smallest distance,
// i.e. distances in non-increasing order.
// [[Rcpp::export]]
IntegerMatrix minimum_spanning_tree(const IntegerVector order) {
  const intx
    n_distances = order.length(),
    n_objects = ceil(sqrt(n_distances + n_distances))
  ;
  unique_ptr<intx[]> left = make_unique<intx[]>(n_distances);
  unique_ptr<intx[]> top = make_unique<intx[]>(n_distances);
Rcout <<"\n\nMST()\n";
  intx k = n_distances;
  for (intx col = n_objects - 1; col--; ) {
    for (intx row = n_objects - 1; row != col; row--) {
      --k;
      left[k] = row;
      top[k] = col;
    }
  }

  unique_ptr<intx[]> island = make_unique<intx[]>(n_objects);
  for (intx i = n_objects; i--; ) {
    island[i] = i;
  }

  IntegerMatrix ret(n_objects - 1, 2);
  intx ret_pos = 0;

  for (intx i = n_distances; i--; ) {
    const intx
      d = order[i],
      left_island = island_housing(left[d], island),
      top_island = island_housing(top[d], island)
    ;
    Rcout << "order[" << i << "] = " << (1+d);
    Rcout << ": " <<(1+ left[d]) << "-"<<(1+top[d]) <<"; ";
    Rcout << "On islands: " << (1+island_housing(left[d], island));
    Rcout << ", " << (1+island_housing(top[d], island)) << ".\n";
    if (top_island != left_island) {
      const intx new_island = (top_island < left_island ? top_island : left_island);
      Rcout << "  - Connecting islands and renumbering to " << new_island <<"\n";
      island[top[d]] = new_island;
      island[left[d]] = new_island;
      island[top_island] = new_island;
      island[left_island] = new_island;
      ret(ret_pos, 0) = top[i] + 1;
      ret(ret_pos, 1) = left[i] + 1;
      ret_pos++;
    }
  }
  return ret;
}
