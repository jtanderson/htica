# Heavy-Tailed ICA Experiments

This directory houses the files used to generate experiments for Heavy-Tailed Independent Component Analysis.
The purpose of the code is to compare various methods of solving the ICA problem when the lower order moments of the source distributions may be undefined.

## Algorithms compared

* [HTICA](http://arxiv.org/abs/1509.00727) - "Heavy-Tailed Independent Component Analysis"
* [FastICA](https://en.wikipedia.org/wiki/FastICA)
  - `pow3` - use fourth cumulant nonlinearity
  - `tanh` - use log cosh nonlinearity
* [Fourier PCA](https://arxiv.org/abs/1306.5825)
* SOBI - _Belouchrani, A., Abed-Meriam, K., Cardoso, J.F. and Moulines, R._ (1997), A blind source separation technique using second-order statistics, IEEE Transactions on Signal Processing, 434–444.
* [JADE](https://en.wikipedia.org/wiki/Joint_Approximation_Diagonalization_of_Eigen-matrices)

## Contents

Following is only a brief description of some key files.
Third party files are located in the `third_party` folder.
See the source for full documentation.

`main.m`: Contains code sections to run each algorithm in sequence, both for dampened and undampened data, then plot the results. This is not a function and is meant to be run directly

`singlecomparison.m`: provides the function `singlecomparison(dim, lowersize, highsize, step, varargin)` which invokes a script to generate heavy-tailed samples, `S`, generates a random mixing matrix `A`, then proceeds to run a chosen ICA algorithm to recover the mixing matrix.
Finally, computes various error metrics and returns them in a 2-row matrix.
Columns of this matrix correspond to each different sample size used for the ICA algorithm

`minkowskiCentroid.m`: uses the Minkowski functional to perform membership queries in the centroid body.
See [1].

`generatesamples.m`: Built to call `mathematicasamples.w` from within Matlab to generate some synethetic data.

`mathematicasamples.w`: A script for the Wolfram kernel to generate samples from nasty Pareto-like distributions.

## References

- 1. J Anderson, N Goyal, A Nandi, L Rademacher. _Heavy-tailed analogues of the covariance matrix for ICA._ Proceedings of the AAAI Conference on Artificial Intelligence 31 (1)
- 2. J Anderson, N Goyal, A Nandi, L Rademacher. _Heavy-tailed independent component analysis_. 2015 IEEE 56th Annual Symposium on Foundations of Computer Science, 290-309
- 3. Belouchrani, A., Abed-Meriam, K., Cardoso, J.F. and Moulines, R. _A blind source separation technique using second-order statistics._ IEEE Transactions on Signal Processing, 434–444.
- 4. Cardoso, Jean-François, and Antoine Souloumiac. _Blind beamforming for non-Gaussian signals._ IEE proceedings F (radar and signal processing). Vol. 140. No. 6. IET Digital Library, 1993.
- 5. Goyal, Navin, Santosh Vempala, and Ying Xiao. _Fourier PCA and robust tensor decomposition._ Proceedings of the forty-sixth annual ACM symposium on Theory of computing. 2014.
- 6. Hyvarinen, A. _Fast ICA for noisy data using Gaussian moments._ 1999 IEEE International Symposium on Circuits and Systems (ISCAS). Orlando, FL, 1999. pp. 57-61 vol.5, doi: 10.1109/ISCAS.1999.777510.

