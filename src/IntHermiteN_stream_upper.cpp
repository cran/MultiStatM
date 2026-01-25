#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List IntHermiteN_stream_upper_cpp(const NumericVector& Fi,
                                        const NumericVector& fi,
                                        const NumericMatrix& Hfi,
                                        int K) {
  int d = Fi.size();
  Rcpp::List res(K);

  long long nrows = 1; // for k = 1: d^0 = 1

  for (int k = 1; k <= K; ++k) {

    if (k == 1) {
      // A[[1]] row is all -1 → CFh = (1 - Fi.x)
      NumericVector out(1);
      double prod = 1.0;
      for (int c = 0; c < d; ++c)
        prod *= (1.0 - Fi[c]);
      out[0] = prod;
      res[k - 1] = out;
    } 
    else {
      nrows *= d;
      NumericVector out(nrows);

      int digits_len = k - 1;
      std::vector<int> digits(digits_len);
      std::vector<int> cnt(d);

      for (long long r = 0; r < nrows; ++r) {

        // Base-d digits of r
        long long tmp = r;
        for (int t = 0; t < digits_len; ++t) {
          digits[t] = tmp % d;
          tmp /= d;
        }

        // Count occurrences
        std::fill(cnt.begin(), cnt.end(), 0);
        for (int t = 0; t < digits_len; ++t)
          cnt[digits[t]]++;

        double prod_cf = 1.0;

        // Build CFh row product
        for (int c = 0; c < d; ++c) {

          int A_val = -1 + cnt[c]; // exactly same A-index logic

          double v;
          if (A_val == -1) {
            v = 1.0 - Fi[c];     // <-- DIFFERENCE FROM LOWER CASE
          } else if (A_val == 0) {
            v = fi[c];
          } else {
            v = Hfi(A_val - 1, c);  // A_val ∈ {1, 2, ..., k-2}
          }

          prod_cf *= v;
        }

        out[r] = prod_cf;
      }

      res[k - 1] = out;
    }
  }

  return res;
}
