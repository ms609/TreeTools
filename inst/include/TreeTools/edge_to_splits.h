#ifndef TreeTools_edge_to_splits_
#define TreeTools_edge_to_splits_

#include <algorithm> /* for sort, unique */
#include <bitset>    /* for count */
#include <climits>   /* for CHAR_BIT */
#include <cstddef>   /* for size_t */
#include <cstdint>   /* for uint64_t */
#include <vector>

namespace TreeTools {

  // Record in `sets` the tips descended from each node, packed
  // `sizeof(Word) * CHAR_BIT` tips per word: node `i` (1-based) occupies
  // words `(i - 1) * n_bin` to `i * n_bin - 1`, and tip `t` is bit
  // `(t - 1) % bits` of word `(t - 1) / bits`.
  // `edge_at(i)` gives the (0-based) index of the `i`th edge to visit; every
  // child must be visited before its parent's own edge.
  // `sets` must hold `n_node * n_bin` zeroed words.
  template <typename Word, typename EdgeAt>
  inline void tips_below(const int* parent, const int* child,
                         const size_t n_edge, const size_t n_tip,
                         const size_t n_bin, EdgeAt edge_at,
                         Word* sets) {
    constexpr size_t bits = sizeof(Word) * CHAR_BIT;
    for (size_t t = 0; t != n_tip; ++t) {
      sets[t * n_bin + t / bits] = Word(1) << (t % bits);
    }
    for (size_t i = 0; i != n_edge; ++i) {
      const size_t e = edge_at(i);
      Word* __restrict__ p = sets + size_t(parent[e] - 1) * n_bin;
      const Word* __restrict__ c = sets + size_t(child[e] - 1) * n_bin;
      for (size_t w = 0; w != n_bin; ++w) {
        p[w] |= c[w];
      }
    }
  }

  // FNV-1a fingerprint of an unrooted topology: a hash of its set of
  // non-trivial splits, so independent of edge order, internal node
  // numbering, and root position.
  // Edges must list each parent before its children, as in preorder or
  // cladewise order; tips are numbered 1 to `n_tip`.
  inline uint64_t topology_hash(const int* parent, const int* child,
                                const size_t n_edge, const size_t n_tip) {
    constexpr uint64_t fnv_offset = 0xcbf29ce484222325ULL;
    constexpr uint64_t fnv_prime = 0x100000001b3ULL;
    const size_t n_bin = (n_tip + 63) / 64;

    size_t n_node = n_tip;
    for (size_t i = 0; i != n_edge; ++i) {
      if (size_t(parent[i]) > n_node) n_node = parent[i];
    }

    std::vector<uint64_t> sets(n_node * n_bin, 0);
    tips_below(parent, child, n_edge, n_tip, n_bin,
               [n_edge](size_t i) { return n_edge - 1 - i; }, sets.data());

    const uint64_t tail_mask = n_tip % 64 ?
      (uint64_t(1) << (n_tip % 64)) - 1 : ~uint64_t(0);

    std::vector<uint64_t> split_hash;
    split_hash.reserve(n_edge);
    for (size_t i = 0; i != n_edge; ++i) {
      if (size_t(child[i]) <= n_tip) continue;
      const uint64_t* split = &sets[size_t(child[i] - 1) * n_bin];
      // Describe each split by the side without tip 1, so that the two
      // halves of an edge hash alike wherever the tree is rooted.
      const uint64_t flip = (split[0] & 1) ? ~uint64_t(0) : 0;
      uint64_t h = fnv_offset;
      size_t n_in = 0;
      for (size_t w = 0; w != n_bin; ++w) {
        uint64_t word = split[w] ^ flip;
        if (w == n_bin - 1) word &= tail_mask;
        n_in += std::bitset<64>(word).count();
        h = (h ^ word) * fnv_prime;
      }
      // Beside a root of degree two, an internal node can carry a trivial
      // split.
      if (n_in > 1 && n_in + 1 < n_tip) split_hash.push_back(h);
    }

    // A split can appear twice when an edge leads to a root of degree two.
    std::sort(split_hash.begin(), split_hash.end());
    split_hash.erase(std::unique(split_hash.begin(), split_hash.end()),
                     split_hash.end());

    uint64_t h = fnv_offset;
    for (const uint64_t s : split_hash) {
      h = (h ^ s) * fnv_prime;
    }
    return h;
  }

}

#endif
