# The Linear Template Fit: *LTF*  

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/badge/version-stable-green.svg)]()
[![](https://img.shields.io/badge/docs-doxygen-blue.svg)](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/) 
[![](https://img.shields.io/badge/docs-home-blue.svg)](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/) 

The *Linear Template Fit* is a matrix formalism for the determination of the best estimator for simulation-based parameter estimation.
The Linear Template Fit combines a linear regression with a linear least-squares method and its optimization. All
predictions need to be calculated beforehand and those are provided for a few selected values for the parameter(s) of interest.
The Linear Template Fit may be useful for a wide variety of fitting problems. More details are provided in [Eur.Phys.J.C 82 (2022) 731](https://doi.org/10.1140/epjc/s10052-022-10581-w) [arXiv:2112.01548](https://arxiv.org/abs/2112.01548).


## Two Implementations
The Linear Template Fit is provided as two (almost) independent implementations, where the difference
is in the linear algebra package that is used for the calculations
- `LTF_Eigen` makes use of the [Eigen](https://eigen.tuxfamily.org) linear algebra package
- `LTF_ROOT` makes use of [ROOT](https://root.cern.ch)'s linear algebra package

Both implementations require a C++17 compatible C++ compiler.

Both implementations provide a set of example programs.


### The Eigen implementaion
The Eigen-based implementation is located in the directory `LTF_Eigen`.
Note that the Eigen-based implementation can be used without `ROOT`. However, for additional functionality, like the generation of plots, this branch can optionally be linked also with ROOT.

### The ROOT implementaion
The ROOT-based implementation is located in the directory `LTF_ROOT`.
The ROOT-based implemenatatoin requires ROOT version 6 (6.26 is recommended) and further requires the tool `CMake3` for compilation.

The ROOT implementation was provided by E. Offermann.


## Links
The pre-print is available from arXiv: [arXiv:2112.01548](https://arxiv.org/abs/2112.01548)

The code repository is hosted at: [github.com/britzger/LinearTemplateFit](https://github.com/britzger/LinearTemplateFit/)

A Doxygen documentation is at: [www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/)

Some lmited further documentation is available from: [www.mpp.mpg.de/~britzger/LinearTemplateFit](https://www.mpp.mpg.de/~britzger/LinearTemplateFit)


## Dependencies
On a machine with CentOS7 and `cvmfs`, one can source `lcgenv-LCG_97-x86_64-centos7-gcc9-opt.sh` to get a recent C++ compiler and ROOT.


