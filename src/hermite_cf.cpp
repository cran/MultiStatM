#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List IntHermiteN_stream_cpp(const NumericVector& Fi,
                                  const NumericVector& fi,
                                  const NumericMatrix& Hfi,
                                  int K) {
  int d = Fi.size();
  Rcpp::List res(K);

  // nrows_k = d^(k-1); we'll keep it incrementally
  long long nrows = 1;  // for k = 1: d^0 = 1

  for (int k = 1; k <= K; ++k) {
    if (k == 1) {
      // A[[1]] all -1 → CFh row = Fi, Si row = all +1
      NumericVector out(1);
      double prodFi = 1.0;
      for (int c = 0; c < d; ++c) {
        prodFi *= Fi[c];
      }
      out[0] = prodFi;
      res[k - 1] = out;
    } else {
      nrows *= d;  // now nrows = d^(k-1)
      NumericVector out(nrows);

      int digits_len = k - 1;
      std::vector<int> digits(digits_len);
      std::vector<int> cnt(d);

      for (long long r = 0; r < nrows; ++r) {
        // 1) base-d digits of r (0-based)
        long long tmp = r;
        for (int t = 0; t < digits_len; ++t) {
          digits[t] = tmp % d;
          tmp /= d;
        }

        // 2) count occurrences of each column index
        std::fill(cnt.begin(), cnt.end(), 0);
        for (int t = 0; t < digits_len; ++t) {
          int j = digits[t];
          if (j >= 0 && j < d)
            cnt[j]++;
        }

        // 3) loop over columns and build products
        double prod_cf = 1.0;
        double prod_si = 1.0;

        for (int c = 0; c < d; ++c) {
          int A_val = -1 + cnt[c];  // A_k[r,c]

          // CFh value
          double v;
          if (A_val == -1) {
            v = Fi[c];
          } else if (A_val == 0) {
            v = fi[c];
          } else {
            // A_val in {1,...,k-2} ⊆ {1,...,K-2}
            v = Hfi(A_val - 1, c);
          }
          prod_cf *= v;

          // Si value: -1 if cnt[c] > 0, else +1
          double s = (cnt[c] > 0) ? -1.0 : 1.0;
          prod_si *= s;
        }

        out[r] = prod_cf * prod_si;
      }

      res[k - 1] = out;
    }
  }

  return res;
}
