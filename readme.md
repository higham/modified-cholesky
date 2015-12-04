`modified-cholesky` - Modified Cholesky factorization
==========

About
-----

`modified-cholesky` contains MATLAB functions that 
compute a modified Cholesky factorization of a symmetric and possibly
indefinite matrix. The algorithm is from

S. H. Cheng and N.J. Higham.
"[A modified Cholesky algorithm based on a symmetric indefinite
  factorization](http://dx.doi.org/10.1137/S0895479896302898)".
SIAM J. Matrix Anal. Appl., 19(4):1097-1110, 1998.

and uses LDL^T factorization with a symmetric form of rook pivoting
proposed by Ashcraft, Grimes, and Lewis.  The functions here are based on
code originally written by Bobby Cheng and Nick Higham in 1996.

The MATLAB functions are:

* `modchol_ldlt`: the modified Cholesky function.  It calls the built-in
  MATLAB function ldl to compute the LDL^T factorization.

* `modchol_ldlt_m`: this is the original version of `modchol_ldlt` from 1996,
  where the `_m` in the name denotes that the LDL^T factorization is computed
  using pure M-code.
  The output of this version should be the same as that from `modchol_ldlt`
  to within rounding error.  The reasons for including this version are as
  follows.
  * Since the code for the factorization is explicitly included as M-code
    the `_m` version is of pedagogical interest.  It will also be useful
    for anyone who wants to modify the factorization to use a different
    pivoting strategy.
  * The `_m` version computes the growth factor for the factorization,
    which this is not available from `modchol_ldlt` itself.

  Note that the `_m` version does not exploit symmetry and is not designed to be
  efficient.  

* `test_modchol_ldlt`: a simple test code.

Requirements
-------------

The codes have been developed under MATLAB 2015b.

License
-------

See `license.txt` for licensing information.
