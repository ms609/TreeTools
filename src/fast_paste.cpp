#include <Rcpp/Lighter>
#include <cstring>
#include <string>
#include <vector>

using namespace Rcpp;

// Format each row of a logical split matrix as "A B | C D".
// `membership` is nSplit x nTip (logical); `labels` is length nTip.
// [[Rcpp::export]]
CharacterVector splits_to_char(const LogicalMatrix membership,
                              const CharacterVector labels) {
  const int n_split = membership.nrow();
  const int n_tip = membership.ncol();
  CharacterVector out(n_split);

  // Cache translated label pointers and lengths
  std::vector<const char *> lab(n_tip);
  std::vector<size_t> lab_len(n_tip);
  for (int j = 0; j < n_tip; ++j) {
    lab[j] = Rf_translateCharUTF8(STRING_ELT(labels, j));
    lab_len[j] = std::strlen(lab[j]);
  }

  std::string buf;
  for (int i = 0; i < n_split; ++i) {
    buf.clear();
    bool first_in = true;
    for (int j = 0; j < n_tip; ++j) {
      if (membership(i, j)) {
        if (!first_in) buf.push_back(' ');
        buf.append(lab[j], lab_len[j]);
        first_in = false;
      }
    }
    buf.append(" | ", 3);
    bool first_out = true;
    for (int j = 0; j < n_tip; ++j) {
      if (!membership(i, j)) {
        if (!first_out) buf.push_back(' ');
        buf.append(lab[j], lab_len[j]);
        first_out = false;
      }
    }
    SET_STRING_ELT(out, i, Rf_mkCharLenCE(buf.c_str(), buf.size(), CE_UTF8));
  }
  return out;
}
