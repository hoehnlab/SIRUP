#include <Rcpp.h>

// Placeholder Rcpp translation unit. The real content (int64-packed DNA
// primitives, PAF parsers, sequence manipulation) lifts over from rad_dev
// in follow-on commits.

//' Confirm the SIRAP shared object loaded
//'
//' Trivial smoke test for the Rcpp build. Returns `TRUE` if the SIRAP
//' shared library is linked and callable.
//'
//' @return `TRUE`
//' @keywords internal
//' @export
// [[Rcpp::export]]
bool sirapRcppLoaded() {
  return true;
}
