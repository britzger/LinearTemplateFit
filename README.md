# The Linear Template Fit: `LTF'  

The *Linear Template Fit* is an analytic expression for parameter estimation which combines the statistical measure with its optimization.
It may be useful for a wide variety of fitting problems. More details are provided in arXiv:2108.XXXXX.

The Linear Template Fit is implemented in `C++` using the [Eigen](https://eigen.tuxfamily.org) package for linear algebra.
Some examples make use of the [ROOT](https://root.cern) analysis framework.
An interactive usability is given through ROOT's [CLING](https://root.cern/cling) LLVM-interpreter for C++ (see below).


If you prefer python, julia, go, awk, or any other language or build-tool, please send me your implementation or the wrapper.
Furthermore, a direct implementation in the ROOT package would be appreciated.


# Links
The pre-print can be downloaded at: arXiv:2108.XXXXX

The code repository is hosted at: [github.com/britzger/LinearTemplateFit](https://github.com/britzger/LinearTemplateFit/)

Doxygen documentation: [www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/)

Possible further documentation: www.mpp.mpg.de/~britzger/LinearTemplateFit](https://www.mpp.mpg.de/~britzger/LinearTemplateFit)



# The LTF package

A doxygen code documentation is at: [mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/)

## Dependencies
The only dependencies are a recent `C++17` compatible C++-compiler and `Eigen` (see below).
The `ROOT`-package is optional.

On CentOS7 with `cvmfs`, you can source `lcgenv-LCG_97-x86_64-centos7-gcc9-opt.sh` to get a recent C++ compiler and ROOT.

## Package structure
The main class of the Linear Template Fit is named `LTF`, and the fit returns a small helper class `LTF::LiTeFit`.
So, the class structure is just:
```
   class LTF
   class LTF::LiTeFit
```

The class `LTF` holds all interface functions and setters for all inputs in different formats.

The class `LTF::LiTeFit` stores all the results as simple members in a standardized format (similar to a C-struct). 
For convenience and historic reasons all (most) members are public, and there is no member naming convention.
All members are correctly instantiatied by `LTF` when calling `LTF::DoLiTeFit()`.

The source code for the shared library `libLTF.so` is 
```
   src/LTF.cxx
   LTF/LTF.h
```


## Eigen
The Linear Template Fit is implemented using the linear algebra package [Eigen](https://eigen.tuxfamily.org).
Eigen is a header-only package and is included into the template fit through the directory `Eigen'.

Please copy a latest release into the directory `LTF/Eigen', or provide a symbolic link.
Alternatively, there is a rather recent version included into this repository.
To make use of it, just generate a symbolic link
```
   cd LTF
   ln -s copy_of_Eigen Eigen
```

## Shared library `libLTF'
Generate a shared library named `libLTF.so` by typing
```
    make slib
````


## Examples and executables

It is recommended to compile the shared library `libLTF.so` before running the examples (to save time and avoid unnecessary recompilation):
```
make slib
```

### Example 1
Example 1 from the writeup is available in two formats: 1) as ROOT-macro, where the distributions are generated with `TRandom3`, or 2) as a standalone executable, where the pseudo-data and the templates are read from the file `data/example_LTF_gaus.txt`.
For details on the Example 1, please read the writeup.

Both examples can be directly compiled, by typing one of the following commands:
```
make example1_LTF_gaus
make example1_LTF_gaus_NoROOT
make all
```
Alternatively, when ROOT is available, the programs can be executed with the ROOT interpreter, using
```
root -b -q -l example1_LTF_gaus.cxx
root -b -q -l example1_LTF_gaus_NoROOT.cxxx
```


### Example 2
The source files of Example 2 are
```
example2_LTF_gaus2D.cxx
example2_LTF_gaus2D_NoROOT.cxx
data/example_LTF_gaus2D.txt
```
Please read Example 1 above for more details on how to compile and execute the example.


### Example 3
Example 3 is only available together with ROOT. Please compile the code or use the ROOT interpreter by calling `make example3_LTF_gaus_sigma` or `root -b -q -l example3_LTF_gaus_sigma.cxx`.


### Example application: strong coupling constant from inclusive jet data
The determination of the strong coupling constant from CMS inclusive jet data is available in the two source files `example_CMSinclusivejets_NN30_BRSSWpaper.cxx` and `example_CMSinclusivejets_MSTW_CMSpaper`.
The required data, covariance matrices, PDF uncrtainties and templates are stored in the directory `data`. 
Execute the example application by compiling the sources or using the ROOT interpreter
```
root -b -q -l example_CMSinclusivejets_MSTW_CMSpaper.cxx
root -b -q -l example_CMSinclusivejets_NN30_BRSSWpaper.cxx
# or:
make example_CMSinclusivejets_NN30_BRSSWpaper
make example_CMSinclusivejets_MSTW_CMSpaper
```

# Plots

A  macro to obtain some plots with ROOT is available in `plot_LTF1D.cxx`.
Though, the plotting routines are certainly not well developed or maintained, but may still serve as an example or for quick studies.
The interface functions are:
```
 void plotLiTeFit(const LTF::LiTeFit& fit, const vector<double>& bins,
                 const string& yaxistitle    = "value [unit]",
                 const string& referencename = "Reference value (#alpha) [unit]",
                 const string& observablename = "Observable [unit]"
```


Include them into the code like
```
#include "plot_LTF1D.cxx"
```
and after the Linear Template Fit pass the object `LTF::LiTeFit` to the macro, like:
```
LTF::LiTeFit fit = ltf.DoLiTeFit();
plotLiTeFit(fit,bins);
```
where `bins` is a `vector<double>` that contains the binning, if this is needed.
The axis titles can be passed to the function as described above.
