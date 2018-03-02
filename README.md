# EPSLASSO

An R package for regression analysis of data from extreme sampling. Both of the low-dimensional (n>p) and high-dimensional (n<p) methods are available.


## Contact
Chao Xu    cxu2@tulane.edu

## Installation
If the R package devtools can be installed, then try:

library("devtools")

install_github("xu1912/EPSLASSO")

If devtools is not available, then have to download the source code and install local using zip (for Windows) or tar.gz (for Linux) file.
Before the installation, check the required packages:

      glmnet(>= 2.0-5),
      
      doParallel (>= 1.0.1),
      
      flare (>= 1.5.0),
      
      foreach (>= 1.4.0),
      
      clime (>= 0.4.1),
      
      mvtnorm (>= 1.0-3),
      
      methods (>= 3.2).
