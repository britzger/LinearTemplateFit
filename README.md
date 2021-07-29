# The Linear Template Fit: `LTF'  

The linear template fit is an analytic expression for the best estimator of a least-square optimization problem
and may become useful for a wide variety of fitting problems. More details are provided in arXiv:2108.XXXXX.

The *Linear Template Fit* is implemented in `C++` using the [Eigen](https://eigen.tuxfamily.org) package for linear algebra.
Some examples make use of the [ROOT](https://root.cern) analysis framework.
An interactive usability is given through ROOT's [CLING](https://root.cern/cling) LLVM-interpreter for C++ (see below).


If you prefer python, julia, go, awk, or any other language or build-tool, please provide us your implementation or the wrapper.

# Links
The pre-print can be downloaded at: arXiv:2108.XXXXX

The code repository is hosted at: [github.com/britzger/LinearTemplateFit](https://github.com/britzger/LinearTemplateFit/)

Doxygen documentation: [www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen](https://www.mpp.mpg.de/~britzger/LinearTemplateFit/doxygen/)

Possible further documentation: www.mpp.mpg.de/~britzger/LinearTemplateFit](https://www.mpp.mpg.de/~britzger/LinearTemplateFit)



# The LTF package

## Package structure
The main class of the Linear Template Fit is named `LTF`, and the fit returns a small helper class `LTF::LiTeFit`.
So, the entire class structure is just:
```
   class LTF
   class LTF::LiTeFit
```

The class `LTF` holds all interface functions and setters for all inputs in different formats.

The class `LTF::LiTeFit` stores all the results as simple members in a standardized format (similar to a C-struct). 
For convenience and historic reasons all (most) members are public, and there is no member naming convention.
All members are correctly instantiatied by `LTF` when calling `LTF::DoLiTeFit()`.

So, the source code for the shared library `libLTF.so` is 
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


