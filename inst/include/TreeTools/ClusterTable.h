#ifndef _TREETOOLS_CLUSTERTABLE_H
#define _TREETOOLS_CLUSTERTABLE_H

#include <bitset> /* for bitset */
#include <vector> /* for vector */
#include "assert.h" /* for ASSERT */
#include "types.h" /* for int16 */
#include "root_tree.h" /* for root_on_node */

#define UNINIT -999
#define INF INTX_MAX

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

#define CT_STACK_SIZE 4

#define CT_IS_LEAF(a) (a) <= n_tip

const int_fast32_t
  CT_MAX_LEAVES = 16383
;

namespace TreeTools {

  class ClusterTable {

    const int16
      L_COL = 0,
      R_COL = 1,
      X_COLS = 2
    ;
    int16
      n_edge,
      n_internal,
      n_leaves,

      n_shared = 0,
      enumeration = 0,
      v_j,
      Tlen,
      Tlen_short,
      Tpos = 0,
      X_ROWS
      ;
    std::vector<int16>
      internal_label,
      leftmost_leaf,
      T,
      visited_nth
      ;
    std::bitset<CT_MAX_LEAVES + 1> Xswitch;
    IntegerMatrix Xarr;

  public:
    ClusterTable(List); // i.e. PREPARE(T)

    inline bool is_leaf(const int16 *v) {
      return *v <= n_leaves;
    }

    inline const int16 edges() {
      return n_edge;
    }

    inline const int16 leaves() {
      return n_leaves;
    }

    inline void ENTER(int16 v, int16 w) {
      T[Tpos++] = v;
      T[Tpos++] = w;
    }

    inline int16 N() {
      return n_leaves;
    }

    inline int16 M() {
      return n_internal;
    }

    inline void TRESET() {
      // This procedure prepares T for an enumeration of its entries,
      // beginning with the first entry.
      Tpos = 0;
    }

    inline void READT(int16 *v, int16 *w) {
      *v = T[Tpos++];
      *w = T[Tpos++];
    }

    inline void NVERTEX(int16 *v, int16 *w) {
      if (Tpos != Tlen) {
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline void NVERTEX_short(int16 *v, int16 *w) {
      // Don't count all-tips or all-ingroup: vertices 0, ROOT, Ingp.
      if (Tpos != Tlen_short) {
        READT(v, w);
        // v_j = *v; // Unneeded unless we go on to call LEFTLEAF
      } else {
        *v = 0;
        *w = 0;
      }
    }

    inline int16 LEFTLEAF() {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return leftmost_leaf[v_j - 1];
    }

    inline void SET_LEFTMOST(int16 index, int16 val) {
      leftmost_leaf[index - 1] = val;
    }

    inline int16 GET_LEFTMOST(int16 index) {
      return leftmost_leaf[index - 1];
    }

    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.

    inline int16 ENCODE(const int16 v) {
      // This function procedure returns as its value the internal label
      // assigned to leaf v
      // MS note: input = v; output = X[v, 3]
      return internal_label[v];
    }

    inline int16 DECODE(const int16 internal_relabeling) {
      // MS: input = X[v, 3], output = v
      return visited_nth[internal_relabeling - 1];
    }

    inline void VISIT_LEAF (const int16* leaf, int16* n_visited) {
      visited_nth[(*n_visited)++] = *leaf;
      internal_label[*leaf] = *n_visited;
    }

    IntegerVector X_decode() {
      IntegerVector ret(N());
      for (int16 i = n_leaves; i--; ) {
        ret(i) = DECODE(i + 1);
      }
      return ret;
    }

    inline int16 X(int16 row, int16 col) {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      return Xarr(col, row - 1);
    }

    inline void setX(int16 row, int16 col, int16 value) {
      ASSERT(row > 0);
      ASSERT(row <= X_ROWS);
      Xarr(col, row - 1) = value;
    }

    IntegerMatrix X_contents() {
      IntegerMatrix ret(X_ROWS, 2);
      for (int16 i = X_ROWS; i--; ) {
        ret(i, 0) = X(i + 1, L_COL);
        ret(i, 1) = X(i + 1, R_COL);
      }
      return ret;
    }

    inline bool CLUSTONL(int16* L, int16* R) {
      return X(*L, L_COL) == *L && X(*L, R_COL) == *R;
    }

    inline bool CLUSTONR(int16* L, int16* R) {
      return X(*R, L_COL) == *L && X(*R, R_COL) == *R;
    }

    inline bool ISCLUST(int16* L, int16* R) {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      return CLUSTONL(L, R) || CLUSTONR(L, R);
    }

    inline void CLEAR() {
      // Each cluster in X has an associated switch that is either cleared or
      // set.
      // This procedure clears every cluster switch in X.
      Xswitch.reset();
    }

    inline void SETSWX(int16* row) {
      Xswitch[*row] = true;
    }

    inline bool GETSWX(int16* row) {
      return Xswitch[*row];
    }

    inline bool NOSWX(const std::size_t& n) {
      return Xswitch.count() == n;
    }

    inline void SETSW(int16* L, int16* R) {
      // If <L,R> is a cluster in X, this procedure sets the cluster switch for <L,R>.
      if (CLUSTONL(L, R)) {
        ++n_shared;
        SETSWX(L);
      } else if (CLUSTONR(L, R)) {
        ++n_shared;
        SETSWX(R);
      }
    }

    inline void UPDATE(){
      // This procedure inspects every cluster switch in X.
      // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R>
      // from X; thereafter ISCLUST(X,L,R) will return the value false.
      for (int16 i = X_ROWS; i--; ) {
        if (!(Xswitch[i])) {
          Xarr(L_COL, i) = 0;
          Xarr(R_COL, i) = 0;
        }
      }
    }

    inline int16 SHARED() {
      return n_shared;
    }

    inline void ADDSHARED() {
      ++n_shared;
    }

    inline void XRESET() {
      // This procedure prepares X for an enumeration of its clusters
      enumeration = 0;
    }

    inline void NCLUS(int16* L, int16* R) {
      // This procedure returns the next cluster <L,R> in the current
      // enumeration of clusters in X.
      // If m clusters are in X, they are returned by the first m invocations
      // of NCLUS after initialization by XRESET; thereafter NCLUS returns the
      // invalid cluster <0,0>.
      *L = X(enumeration, 0);
      *R = X(enumeration, 1);
      ++enumeration;
    }

  };

  ClusterTable::ClusterTable(List phylo) {

    const List rooted = TreeTools::root_on_node(phylo, 1);
    const IntegerMatrix edge = rooted["edge"];

    // BEGIN
    n_internal = rooted["Nnode"]; // = M
    CharacterVector leaf_labels = rooted["tip.label"];
    if (leaf_labels.length() > CT_MAX_LEAVES) {
      throw std::length_error("Tree has too many leaves. "
                                "Contact the 'TreeTools' maintainer.");
    }
    n_leaves = leaf_labels.length(); // = N
    n_edge = edge.nrow();
    const int16 n_vertex = M() + N();
    Tlen = 2 * n_vertex;
    Tlen_short = Tlen - (2 * 3);
    T = std::vector<int16> (Tlen);

    leftmost_leaf = std::vector<int16> (n_vertex);
    visited_nth = std::vector<int16> (n_leaves);
    internal_label = std::vector<int16>(1 + n_leaves); // We're not using -1.
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
      parent_i = edge(i, 0),
        child_i = edge(i, 1)
      ;
      if (!GET_LEFTMOST(parent_i)) {
        SET_LEFTMOST(parent_i, GET_LEFTMOST(child_i));
      }
      if (is_leaf(&child_i)) {
        VISIT_LEAF(&child_i, &n_visited);
        ++weights[parent_i];
        ENTER(child_i, 0);
      } else {
        weights[parent_i] += 1 + weights[child_i];
        ENTER(child_i, weights[child_i]);
      }
    }
    ENTER(edge(0, 0), weights[edge(0, 0)]);

    // BUILD Cluster table
    X_ROWS = n_leaves;
    Xarr = IntegerMatrix(X_COLS, X_ROWS);
    // Xswitch = std::bitset<DAY_MAX_LEAVES>;

    // This procedure constructs in X descriptions of the clusters in a
    // rooted tree described by the postorder sequence T with weights,
    // BUILD assigns each leaf an internal label so that every cluster
    // is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
    // simply described by a pair <L,R> of internal labels.

    TRESET();
    for (int16 i = 1; i != N(); ++i) {
      setX(i, L_COL, 0);
      setX(i, R_COL, 0);
    }
    int16 leafcode = 0, v, w, L, R = UNINIT, loc;

    NVERTEX(&v, &w);
    while (v) {
      if (is_leaf(&v)) {
        ++leafcode;
        // We prepared the encoder in an earlier step, so need no X[v, 3] <- leafcode
        R = leafcode;
        NVERTEX(&v, &w);
      } else {
        L = ENCODE(LEFTLEAF());
        NVERTEX(&v, &w);
        loc = w == 0 ? R : L;
        setX(loc, L_COL, L);
        setX(loc, R_COL, R);
      }
    }
  }
}

#endif
