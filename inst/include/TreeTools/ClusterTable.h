#ifndef _TREETOOLS_CLUSTERTABLE_H
#define _TREETOOLS_CLUSTERTABLE_H

#include <array> /* for array */
#include <bitset> /* for bitset */
#include <vector> /* for vector */
#include <Rcpp/Lightest>
#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16 */
#include "root_tree.h" /* for root_on_node */

#define UNINIT -999
#define INF TreeTools::INTX_MAX

#define CT_PUSH(a, b, c, d)                                      \
  S[Spos++] = (a);                                               \
  S[Spos++] = (b);                                               \
  S[Spos++] = (c);                                               \
  S[Spos++] = (d)

#define CT_POP(a, b, c, d)                                       \
  (d) = S[--Spos];                                               \
  (c) = S[--Spos];                                               \
  (b) = S[--Spos];                                               \
  (a) = S[--Spos]

// Required by TreeDist 2.9.2
// TODO Remove in later version, to prefer ct_stack_size
#define CT_STACK_SIZE 4

#define CT_IS_LEAF(a) (a) <= n_tip

// Required by TreeDist 2.9.2
// TODO Remove in later version, to prefer ct_max_leaves
const int_fast32_t CT_MAX_LEAVES = 16383;

namespace TreeTools {

  constexpr int_fast32_t ct_max_leaves = 16383;
  constexpr int_fast32_t ct_stack_size = 4;
  
  
  class ClusterTable {
    
    struct ClusterRow {
      int16 L;
      int16 R;
    };
    
    int16 n_edge;
    int16 n_internal;
    int16 n_leaves;
    int16 n_shared = 0;
    int16 enumeration = 0;
    int16 v_j;
    int16 Tlen;
    int16 Tlen_short;
    int16 X_ROWS;
    std::vector<int16> internal_label;
    int16 *internal_label_ptr = nullptr;
    std::vector<int16> leftmost_leaf;
    std::vector<int16> T;
    int16* T_ptr = nullptr;
    std::vector<int16> visited_nth;
    std::vector<ClusterRow> x_rows;
    // Using bitset; can obtain a ~1% speedup using vector of ULLs
    // Retaining slower code as easier to read.
    // See branch ct-xswitch for implementation
    std::bitset<CT_MAX_LEAVES + 1> Xswitch;
    // Track number of set switches (excluding index 0)
    std::size_t xswitch_set_count = 0;
    

  public:
    ClusterTable(Rcpp::List); // i.e. PREPARE(T)

    [[nodiscard]] inline bool is_leaf(const int16 v) noexcept {
      return v <= n_leaves;
    }
    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer is_leaf(int16 v)
    [[nodiscard]] inline bool is_leaf(const int16 *v) noexcept {
      return *v <= n_leaves;
    }

    [[nodiscard]] inline const int16 edges() noexcept {
      return n_edge;
    }

    [[nodiscard]] inline const int16 leaves() noexcept {
      return n_leaves;
    }

    inline void ENTER(int16 v, int16 w) noexcept {
      *T_ptr = v;
      ++T_ptr;
      *T_ptr = w;
      ++T_ptr;
    }

    [[nodiscard]] inline int16 N() noexcept {
      return n_leaves;
    }

    [[nodiscard]] inline int16 M() noexcept {
      return n_internal;
    }

    inline void TRESET() noexcept {
      // This procedure prepares T for an enumeration of its entries,
      // beginning with the first entry.
      T_ptr = T.data();
    }

    inline void READT(int16 *v, int16 *w) {
      *v = *T_ptr++;
      *w = *T_ptr++;
    }

    inline void NVERTEX(int16 *v, int16 *w) noexcept {
      if (T_ptr != T.data() + Tlen) {
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline void NVERTEX_short(int16 *v, int16 *w) noexcept {
      // Don't count all-tips or all-ingroup: vertices 0, ROOT, Ingp.
      if (T_ptr != T.data() + Tlen_short) {
        READT(v, w);
        // v_j = *v; // Unneeded unless we go on to call LEFTLEAF
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline int16 LEFTLEAF() noexcept {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return leftmost_leaf[v_j - 1];
    }

    inline void SET_LEFTMOST(int16 index, int16 val) noexcept {
      leftmost_leaf[index - 1] = val;
    }

    [[nodiscard]] inline int16 GET_LEFTMOST(int16 index) noexcept {
      return leftmost_leaf[index - 1];
    }

    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.

    inline int16 ENCODE(const int v) noexcept {
      // This function procedure returns as its value the internal label
      // assigned to leaf v
      // MS note: input = v; output = X[v, 3]
      return internal_label_ptr[v];
    }

    inline int16 DECODE(const int16 internal_relabeling) noexcept {
      // MS: input = X[v, 3], output = v
      return visited_nth[internal_relabeling - 1];
    }

    inline void VISIT_LEAF (const int16* leaf, int16* n_visited) {
      visited_nth[(*n_visited)++] = *leaf;
      internal_label[*leaf] = *n_visited;
    }

    Rcpp::IntegerVector X_decode() {
      Rcpp::IntegerVector ret(N());
      for (int16 i = n_leaves; i--; ) {
        ret(i) = DECODE(i + 1);
      }
      return ret;
    }
    
    [[nodiscard]] inline int16 X_left(int16 row) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      ASSERT(x_rows[row - 1].L < std::numeric_limits<int16>::max());
      return x_rows[row - 1].L;
    }
    
    [[nodiscard]] inline int16 X_right(int16 row) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      ASSERT(x_rows[row - 1].R < std::numeric_limits<int16>::max());
      return x_rows[row - 1].R;
    }
    
    inline void setX_left(int16 row, int16 value) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      x_rows[row - 1].L = value;
    }
    
    inline void setX_right(int16 row, int16 value) noexcept {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      x_rows[row - 1].R = value;
    }
    
    Rcpp::IntegerMatrix X_contents() noexcept {
      Rcpp::IntegerMatrix ret(X_ROWS, 2);
      for (int16 i = X_ROWS; i--; ) {
        ret(i, 0) = x_rows[i].L;
        ret(i, 1) = x_rows[i].R;
      }
      return ret;
    }

    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer CLUSTONL(int16 L, R)
    [[nodiscard]] inline bool CLUSTONL(int16* L, int16* R) noexcept {
      return X_left(*L) == *L && X_right(*L) == *R;
    }
    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer CLUSTONR(int16 L, R)
    [[nodiscard]] inline bool CLUSTONR(int16* L, int16* R) noexcept {
      return X_left(*R) == *L && X_right(*R) == *R;
    }
    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer ISCLUST(int16 L, R)
    [[nodiscard]] inline bool ISCLUST(int16* L, int16* R) noexcept {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      return CLUSTONL(*L, *R) || CLUSTONR(*L, *R);
    }

    [[nodiscard]] inline bool CLUSTONL(int16 L, int16 R) noexcept {
      ASSERT(L > 0 && L <= X_ROWS);
      const ClusterRow &r = x_rows[L - 1];
      return (r.L == L) & (r.R == R);
    }
    
    [[nodiscard]] inline bool CLUSTONR(int16 L, int16 R) noexcept {
      ASSERT(R > 0 && R <= X_ROWS);
      const ClusterRow &r = x_rows[R - 1];
      return (r.L == L) & (r.R == R);
    }
    
    [[nodiscard]] inline bool ISCLUST(int16 L, int16 R) noexcept {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      ASSERT(L > 0 && L <= X_ROWS);
      const ClusterRow &r_L = x_rows[L - 1];
      if ((r_L.L == L) & (r_L.R == R)) return true;
      
      ASSERT(L != R);
      ASSERT(R > 0 && R <= X_ROWS);
      const ClusterRow &r_R = x_rows[R - 1];
      return (r_R.L == L) & (r_R.R == R);
    }

    inline void CLEAR() noexcept {
      // Each cluster in X has an associated switch that is either cleared or
      // set.
      // This procedure clears every cluster switch in X.
      Xswitch.reset();
      xswitch_set_count = 0;
    }
    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer SETSWX(int16 row)
    inline void SETSWX(int16* row) noexcept {
      // Only increment our counter on a 0 -> 1 transition
      const auto idx = static_cast<std::size_t>(*row);
      if (!Xswitch[idx]) {
        Xswitch[idx] = true;
        ++xswitch_set_count;
      }
    }
    
    inline void SETSWX(std::size_t row) noexcept {
      // Only increment our counter on a 0 -> 1 transition
      if (!Xswitch[row]) {
        Xswitch[row] = true;
        ++xswitch_set_count;
      }
    }

    [[nodiscard]] inline bool GETSWX(int16* row) noexcept {
      return Xswitch[*row];
    }

    [[nodiscard]] inline bool NOSWX(const std::size_t& n) noexcept {
      return xswitch_set_count == n;
    }
    
    // Required by TreeDist 2.9.2
    // TODO Remove in later version, to prefer SETSW(int16 L, R)
    inline void SETSW(int16 *L, int16 *R) noexcept {
      const int16 l = *L;
      const int16 r = *R;
      // If <L,R> is a cluster in X, 
      // this procedure sets the cluster switch for <L,R>.
      if (CLUSTONL(l, r)) {
        ++n_shared;
        SETSWX(l);
      } else if (CLUSTONR(l, r)) {
        ++n_shared;
        SETSWX(r);
      }
    }
    
    inline void SETSW(int16 L, int16 R) noexcept {
      // If <L,R> is a cluster in X, 
      // this procedure sets the cluster switch for <L,R>.
      if (CLUSTONL(L, R)) {
        ++n_shared;
        SETSWX(L);
      } else if (CLUSTONR(L, R)) {
        ++n_shared;
        SETSWX(R);
      }
    }

    inline void UPDATE() noexcept {
      // This procedure inspects every cluster switch in X.
      // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R>
      // from X; thereafter ISCLUST(X,L,R) will return the value false.
      for (int16 i = X_ROWS; i--; ) {
        if (!(Xswitch[i])) {
          x_rows[i].L = 0;
          x_rows[i].R = 0;
        }
      }
    }

    [[nodiscard]] inline int16 SHARED() noexcept {
      // Used by COMCLUST in TreeDist::Day_1985.cpp
      return n_shared;
    }

    inline void ADDSHARED() noexcept {
      ++n_shared;
    }

    inline void XRESET() noexcept {
      // This procedure prepares X for an enumeration of its clusters
      enumeration = 0;
    }

    inline void NCLUS(int16* L, int16* R) noexcept {
      // This procedure returns the next cluster <L,R> in the current
      // enumeration of clusters in X.
      // If m clusters are in X, they are returned by the first m invocations
      // of NCLUS after initialization by XRESET; thereafter NCLUS returns the
      // invalid cluster <0,0>.
      ASSERT(enumeration < X_ROWS);
      const int16 row = static_cast<int16>(enumeration + 1); // rows are 1-based
      *L = X_left(row);
      *R = X_right(row);
      ++enumeration;
    }

  };

  inline ClusterTable::ClusterTable(Rcpp::List phylo) {

    const Rcpp::List rooted = TreeTools::root_on_node(phylo, 1);
    const Rcpp::IntegerMatrix edge = rooted["edge"];

    // BEGIN
    n_internal = rooted["Nnode"]; // = M
    Rcpp::CharacterVector leaf_labels = rooted["tip.label"];
    if (leaf_labels.length() > int(CT_MAX_LEAVES)) {
      Rcpp::stop("Tree has too many leaves. "
                 "Contact the 'TreeTools' maintainer.");
    }
    ASSERT(CT_MAX_LEAVES <= std::numeric_limits<int16>::max());
    n_leaves = int16(leaf_labels.length()); // = N
    if (double(edge.nrow()) > double(std::numeric_limits<int16>::max())) {
      Rcpp::stop("Tree has too many edges. "
                 "Contact the 'TreeTools' maintainer.");
    }
    n_edge = int16(edge.nrow());
    const int16 n_vertex = M() + N();
    Tlen = 2 * n_vertex;
    Tlen_short = Tlen - (2 * 3);
    T = std::vector<int16>(Tlen);
    T_ptr = T.data();

    leftmost_leaf.reserve(n_vertex);
    visited_nth.reserve(n_leaves);
    internal_label.reserve(1 + n_leaves); // We're not using -1.
    internal_label_ptr = internal_label.data();
    int16 n_visited = 0;
    std::vector<int16> weights(1 + n_vertex);

    for (int16 i = 1; i != n_leaves + 1; ++i) {
      SET_LEFTMOST(i, i);
      weights[i] = 0;
    }
    for (int16 i = 1 + n_leaves; i != 1 + n_vertex; ++i) {
      SET_LEFTMOST(i, 0);
      weights[i] = 0;
    }
    for (int16 i = n_edge; i--; ) {
      const int16
        parent_i = int16(edge(i, 0)),
        child_i = int16(edge(i, 1))
      ;
      if (!GET_LEFTMOST(parent_i)) {
        SET_LEFTMOST(parent_i, GET_LEFTMOST(child_i));
      }
      if (is_leaf(child_i)) {
        VISIT_LEAF(&child_i, &n_visited);
        ++weights[parent_i];
        ENTER(child_i, 0);
      } else {
        weights[parent_i] += 1 + weights[child_i];
        ENTER(child_i, weights[child_i]);
      }
    }
    ENTER(int16(edge(0, 0)), weights[edge(0, 0)]);

    // BUILD Cluster table
    X_ROWS = n_leaves;
    x_rows = std::vector<ClusterRow>(X_ROWS);

    // This procedure constructs in X descriptions of the clusters in a
    // rooted tree described by the postorder sequence T with weights,
    // BUILD assigns each leaf an internal label so that every cluster
    // is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
    // simply described by a pair <L,R> of internal labels.

    TRESET();
    for (int16 i = 1; i != N(); ++i) {
      setX_left(i, 0);
      setX_right(i, 0);
    }
    int16 leafcode = 0, v, w, L, R = UNINIT, loc;

    NVERTEX(&v, &w);
    while (v) {
      if (is_leaf(v)) {
        ++leafcode;
        // We prepared the encoder in an earlier step, so need no X[v, 3] <- leafcode
        R = leafcode;
        NVERTEX(&v, &w);
      } else {
        L = ENCODE(LEFTLEAF());
        NVERTEX(&v, &w);
        loc = w == 0 ? R : L;
        setX_left(loc, L);
        setX_right(loc, R);
      }
    }
  }
}

#endif
