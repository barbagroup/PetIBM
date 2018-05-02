# Use PetIBM API

In the installation path, PetIBM also provides a library and header files that 
have all necessary API to build a flow solver.

## API Documentation

The API documentation can be found in the **Modules** section of the 
[Doxygen pages](https://barbagroup.github.io/PetIBM/modules.html). The subsection
of **PetIBM building blocks** shows the details of the elements for building
a flow solver. And the subsection of **Flow solvers and utilities** shows the
flow solvers built on those building blocks.

A local version of Doxygen pages can be created by `doxygen` utility. Go into
the **doc** folder of PetIBM, and run `doxygen` command. HTML pages will be
created and put into a subfolder called **html**. Users then can open up the 
main page of the manual, **html/index.html**, using a preferred browser. 

## API usage examples

The source code of the flow solvers mentioned in the API documentation may be
unfriendly for beginners, so we provide other simplified solvers in the folder
**$PETIBM_SOURCE/examples/api_examples**. 
The one in the sub-folder **liddrivencavity2d**
provides step-by-step instruction in the comments of the source code. See the
README files in those folders.

## Compiling and linking against PetIBM

### Header files

Let **$PETIBM_INSTALL** be the installation path of PetIBM. All header files for
public PetIBM API are under **$PETIBM_INSTALL/include/petibm**, while those for
third-party dependencies, if any, are under **$PETIBM_INSTALL/include**. For example,
if users use `--enable-amgxwrapper` flag while configuring PetIBM, the 
`AmgXSolver.hpp` will show up in **$PETIBM_INSTALL/include**.

There are many headers of PetIBM, while only the top-level ones are normally
used in an application code:

* `type.h`: useful data types in PetIBM
* `parser.h`: functions to parse input YAML file
* `mesh.h`: abstract mesh object and its factory function
* `boundary.h`: abstract boundary (and ghost points) object and its factory function
* `solution.h`: abstract solution object and its factory function 
* `operators.h`: factory functions for creating operators (matrices)
* `linsolver.h`: abstract object for linear solvers and its factory function
* `bodypack.h`: abstract object for multi-body holders and its factory function
* `timeintegration.h`: abstract object for temporal integration and its factory function

To compile with PetIBM headers, remember to have `-I$PETIBM_INSTALL/include` in the
compilation flags. Also, the path to PETSc headers is necessary, e.g.,
`-I$PETSC_DIR/PETSC_ARCH/include`.

If users provide their own installations of third-party dependencies to build 
PetIBM, remember to also include those paths. For example, if users provide 
`--with-yamlcpp-dir=$YAMLCPP_DIR` during configuration, then they will have to 
add `-I$YAMLCPP_DIR/include` to the compilation flags.

### Library

PetIBM provides a single library, either **libpetibm.a** (static library) or 
**libpetibm.so** (shared library) to link against. The library is under
**$PETIBM_INSTALL/lib**. Use `-L$PETIBM_INSTALL/lib -lpetibm` to link.

There may be other libraries in **$PETIBM_INSTALL/lib** as dependencies
of PetIBM, if users let PetIBM download and build them during configuration.
For example, if users use `--enable-yamlcpp` during configuring PetIBM, the 
yaml-cpp library will show up. And hence users can link against those libraries
from **$PETIBM_INSTALL/lib**. Note, in some systems or OS, the third-party 
libraries may be in **$PETIBM_INSTALL/lib64**.

Also, linking against PETSc library is necessary. The flag for linking against
PETSc can be found in the file **$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables**.
Look for the line starting with `PETSC_WITH_EXTERNAL_LIB = `.

An alternative way to link PetIBM and let PetIBM library automatically handle those
dependencies is through `libtool`. PetIBM library is built using `libtool`, and
`libtool` records the information of dependencies in an abstract layer
**$PETIBM_INSTALL/lib/libpetibm.la**. Users can see the information of dependencies
in the `.la` file, but the proper way to use `libtool` and the `.la` file is

```
libtool --tag=CXX --mode=link $CXX -o <executable name> <object files> \
    $LINKER_FLAGS $PETIBM_INSTALL/lib/libpetibm.la
```

Or if doing compiling and linking at once

```
libtool --tag=CXX --mode=link $CXX -o <executable name> <source files>  \
    $CXX_FLAGS $LINKER_FLAGS $PATH_TO_HEADERS $PETIBM_INSTALL/lib/libpetibm.la
```

And then `libtool` will handle the dependencies.
