#include <Rcpp/Lightest>
#include <vector> /* for vector */
#include <cmath> /* for ceil, sqrt() */
#include "../inst/include/TreeTools/types.h"
using namespace Rcpp;
using namespace std;

int16 island_housing(int16 x, vector<int16> &island) {
  if (x == island[x]) return x;
  int16 parent_house = island_housing(island[x], island);
  if (parent_house == island[x]) return island[x];
  island[x] = parent_house;
  return parent_house;
}

// @param order Index of largest, second largest, ... smallest distance,
// i.e. distances in non-increasing order.
// [[Rcpp::export]]
IntegerMatrix minimum_spanning_tree(const IntegerVector order) {
  const int32 n_distances = order.length();
  const int16 n_objects = ceil(sqrt(n_distances + n_distances));
  vector<int16> left (n_distances);
  vector<int16> top (n_distances);

  int32 k = n_distances;
  for (int16 col = n_objects - 1; col--; ) {
    for (int16 row = n_objects - 1; row != col; row--) {
      --k;
      left[k] = row;
      top[k] = col;
    }
  }

  vector<int16> island (n_objects);
  for (int16 i = n_objects; i--; ) {
    island[i] = i;
  }

  IntegerMatrix ret(n_objects - 1, 2);
  int16 ret_pos = 0;

  for (int32 i = n_distances; i--; ) {
    const int32 d = order[i];
    const int16
      left_island = island_housing(left[d], island),
      top_island = island_housing(top[d], island)
    ;
    if (top_island != left_island) {
      const int16 new_island = (top_island < left_island ? top_island : left_island);
      island[top[d]] = new_island;
      island[left[d]] = new_island;
      island[top_island] = new_island;
      island[left_island] = new_island;
      ret(ret_pos, 0) = top[d] + 1;
      ret(ret_pos, 1) = left[d] + 1;
      if (ret_pos == n_objects - 2) break;
      ret_pos++;
    }
  }
  return ret;
}
