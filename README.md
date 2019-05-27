2x2 SVD of a matrix (Singular Value Decomposition) (C++) inspired from Eigen.

(MPL 2.0-license)

<link href="https://fonts.googleapis.com/css?family=Cookie" rel="stylesheet"><a class="bmc-button" target="_blank" href="https://www.buymeacoffee.com/jBnA3c2Fw"><img src="https://www.buymeacoffee.com/assets/img/BMC-btn-logo.svg" alt="Buy me a coffee"><span style="margin-left:5px">You can buy me a coffee</span></a> if you want to support o(^◇^)o.

Usage

```
     #include "svd_2x2.hpp"
     Tbx::Mat2 M( a, b,
                  c, d);
     // Compute full svd of M 
     Tbx::SVD_2x2 svd(M);

     // Retrieve M = U * S * Vt
     Tbx::Mat2 Vt = svd.matrix_v().transpose();
     Tbx::Mat2 U = svd.matrix_u();
     Tbx::Vec2 = svd.singular_values(); // Singular values diagonal matrix as a simple vector.
```
