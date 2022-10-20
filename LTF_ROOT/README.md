# The Linear Template Fit: *LTF*  

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/badge/version-stable-green.svg)]()
[![](https://img.shields.io/badge/docs-doxygen-blue.svg)](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/) 
[![](https://img.shields.io/badge/docs-home-blue.svg)](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/) 

The *Linear Template Fit* is a matrix formalism for the determination of the best estimator for simulation-based parameter estimation.
The Linear Template Fit combines a linear regression with a linear least-squares method and its optimization. All
predictions need to be calculated beforehand and those are provided for a few selected values for the parameter(s) of interest.
The Linear Template Fit may be useful for a wide variety of fitting problems. More details are provided in [Eur.Phys.J.C 82 (2022) 731](https://doi.org/10.1140/epjc/s10052-022-10581-w) [arXiv:2112.01548](https://arxiv.org/abs/2112.01548).


## ROOT Implementation
The present implementation makes use of the [linear algebra classes](https://root.cern.ch/root/htmldoc/guides/users-guide/LinearAlgebra.html) as implemented in the [ROOT](root.cern.ch) analysis framework.
Hence, the only external dependency is ROOT.

Several examples are available in the example diretory and are compiled together with the main library.
Besides the numerical results, the examples provide also some plots based on the ROOT plotting tools.

All examples and classes can directly be executed in ROOT's C++-interpreter [CLING](https://root.cern/cling) (see instructions below).

If you use this code please cite: *D. Britzger, "The Linear Template Fit", [Eur.Phys.J.C 82 (2022) 731](https://doi.org/10.1140/epjc/s10052-022-10581-w) [[arXiv:2112.01548]](https://arxiv.org/abs/2112.01548), [DOI:10.1140/epjc/s10052-022-10581-w](https://doi.org/10.1140/epjc/s10052-022-10581-w)*.

The ROOT implementation was provided by E. Offermann.


## Compilation
The compilation requires a C++17 compatible compiler and CMake version 3.
In order to compile the package and the examples, type in the main directory (not in the directory `LTF`):
```
cmake -S . -B build
cmake --build build
```

After successful compilation, you see following files:
```
build/lib:
libLTF.rootmap
libLTF_rdict.pcm
libLTF.so

build/bin:
example1_LTF_gaus
example2_LTF_gaus2D
example_CMSinclusivejets_MSTW_CMSpaper
example3_LTF_gaus_sigma
example_CMSinclusivejets_NN30_BRSSWpaper
```


Now, you can directly execute the examples:
```
build/bin/example1_LTF_gaus
build/bin/example2_LTF_gaus2D
build/bin/example3_LTF_gaus_sigma

```
Some of these examples print several plots into a new directory, called `plots`.


The additional examples of the alpha_s fit from the paper require the directory `data` to be directly accessible.
Thus, first call once:
```
ln -s examples/data .
```
Then you can directly execute the alpha_s-fits from the (LTF paper)[https://doi.org/10.1140/epjc/s10052-022-10581-w]
```        
build/bin/example_CMSinclusivejets_MSTW_CMSpaper
build/bin/example_CMSinclusivejets_MSTW_CMSpaper
```

## Interactive usabilitiy
Root's interactive LLV-interpreter `CLING` can directly interpret the source code without explicitly calling a constructor.
Just call
```
$ root -b
root [0] .I ./LTF/inc
root [1] .L ./LTF/src/LTF.cxx 
root [2] .L ./LTF/src/LTF_Tools.cxx 
root [3] .L ./LTF/src/LTF_ROOTTools.cxx
root [4] .x examples/example1_LTF_gaus.cxx
```

Another option to directly load the shared library would be
```
$ root 
root [0] gSystem->AddIncludePath("-I./build/include")        // or: .I build/include
root [1] #include "LTF/LTF.h"
root [2] gSystem->Load("build/lib/libLTF.so")
root [3] LTF ltf
// continue to work...
```

Alternative options would be, to modify the `.rootlogon.C` (commonly in `$HOME`) an/or add
```
 export LD_LIBRARY_PATH=$PWD/build/lib:$LD_LIBRARY_PATH
 export PATH=$PWD/build/bin:$LD_LIBRARY_PATH
```


## Further Links
The pre-print is available from arXiv: [arXiv:2112.01548](https://arxiv.org/abs/2112.01548)

The code repository is hosted at: [github.com/britzger/LinearTemplateFit](https://github.com/britzger/LinearTemplateFit/)

A Doxygen documentation is at: [www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/)

Some limited further documentation is available from: [www.mpp.mpg.de/~britzger/LinearTemplateFit](https://www.mpp.mpg.de/~britzger/LinearTemplateFit)

The ROOT analysis framework: [root.cern.ch](root.cern.ch)


## Package structure
The main class of the Linear Template Fit is named `LTF`, and the fit returns a small helper class `LTF::LiTeFit`.
So, the class structure is simply:
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
   inc/LTF/LTF.h
```


## Helper functions
Some helper functions are available to read files or for basic plotting:
```
LTF_Tools.h
LTF_ROOTTools.h
plot_LTF1D.cxx
```


