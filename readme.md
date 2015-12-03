`modified-cholesky` - Modified Cholesky factorization
==========

About
-----

`modified-cholesky` contains the MATLAB function `modchol_ldlt`, which
computes a modified Cholesky factorization of a symmetric and possibly
indefinite matrix. The algorithm is from

S. H. Cheng and N.J. Higham.
"[A modified Cholesky algorithm based on a symmetric indefinite
  factorization](http://dx.doi.org/10.1137/S0895479896302898)".
SIAM J. Matrix Anal. Appl., 19(4):1097-1110, 1998.

and it uses LDL^T factorization with a symmetric form of rook pivoting
proposed by Ashcraft, Grimes, and Lewis.

The MATLAB functions are:

* `modchol_ldlt``: the modified Cholesky function.  
This is a research code that does not exploit symmetry and is not
designed to be efficient.

* `test_modchol_ldlt``: a simple test code.

Requirements
-------------

The codes have been developed under MATLAB 2015b.

License
-------

See `license.txt` for licensing information.
