# LinearTemplateFit
The Linear Template Fit



# Install
The *Linear Template Fit* is implemented in `C++` using the [Eigen](https://eigen.tuxfamily.org) package for linear algebra.
Some examples make use of the [ROOT](https://root.cern) analysis framework.

# Package structure
The main class of the Linear Template Fit is named `LTF`, and the fit returns a small helper class `LTF::LiTeFit`.
So, the entire class structure is just:
```
   class LTF
   class LTF::LiTeFit
```

The `class LTF` keeps all interface functions and setters for different input formats.

The `class LTF::LiTeFit` keeps all the results as simple members in a standardized format (like a `C`-struct). 
For convenience and historic reasons all (most) members are public without any naming convention, but 
all members are correctly instantiatied by `LTF` when calling `LTF::DoLiTeFit()`.


## Eigen
The Linear Template Fit is implemented using the linear algebra package [Eigen](https://eigen.tuxfamily.org).
Eigen is a header-only package and is included into the template fit through the directory `Eigen'.

Please copy a latest release into the directory `LTF/Eigen', or provide a symbolic link.
Alternatively, there is a rather recent stub version included into this repository.
To make use of it, just generate a symbolic link
```
   cd LTF
   ln -s copy_of_Eigen Eigen
```

## Shared library libLTF
Generate a shared library named `libLTF.so` by typing
```
    make slib
````

