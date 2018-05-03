# Use PetIBM API

In the installation path, PetIBM also provides a library and header files that have all necessary API to build a flow solver.

## API Documentation

The API documentation can be found in the **Modules** section of the [Doxygen pages](https://barbagroup.github.io/PetIBM/modules.html).
The subsection of **PetIBM building blocks** shows details of the elements for building a flow solver.
The subsection of **Flow solvers and utilities** shows the flow solvers built on those building blocks.

A local version of Doxygen pages can be created with the `doxygen` utility.
To do that, go into the `doc` folder of PetIBM, and run the `doxygen` command.
HTML pages will be created in the subfolder `html` of the `doc` directory.
Users can, then, open the main page of the manual, `html/index.html`, using their preferred browser. 

## API usage examples

The source code of the flow solvers mentioned in the API documentation may be unfriendly for beginners, so we provide other simplified solvers in the folder `$PETIBM_DIR/examples/api_examples`.
For example, `liddrivencavity2d` provides step-by-step instructions in the comments of the source code.
(See the README files in those folders for additional details.)

## Compiling and linking against PetIBM

### Header files

Let `$PETIBM_INSTALL` be the installation path of PetIBM.
All header files for public PetIBM API are under `$PETIBM_INSTALL/include/petibm`, while those for third-party dependencies, if any, are under `$PETIBM_INSTALL/include`.
For example, if users configure PetIBM using the flag `--enable-amgxwrapper`, the header file `AmgXSolver.hpp` will be copied to `$PETIBM_INSTALL/include`.

There are many headers of PetIBM, while only the top-level ones are normally
used in an application code:

* `type.h`: useful data types in PetIBM,
* `parser.h`: functions to parse YAML input file,
* `mesh.h`: abstract mesh object and its factory functions,
* `boundary.h`: abstract boundary (and ghost points) object and its factory functions,
* `solution.h`: abstract solution object and its factory functions,
* `operators.h`: factory functions for creating operators (matrices),
* `linsolver.h`: abstract object for linear solvers and its factory functions,
* `bodypack.h`: abstract object for multi-body holders and its factory functions,
* `timeintegration.h`: abstract object for temporal integration and its factory functions.

To compile with PetIBM headers, remember to have `-I$PETIBM_INSTALL/include` in the compilation flags.
Also, the path to PETSc headers is necessary, e.g., `-I$PETSC_DIR/PETSC_ARCH/include`.

If users provide their own installations of third-party dependencies to build 
PetIBM, remember to also include those paths.
For example, if users provide `--with-yamlcpp-dir=$YAMLCPP_DIR` during configuration, then they will have to add `-I$YAMLCPP_DIR/include` to the compilation flags.

### Library

PetIBM provides a single library, either `libpetibm.a` (static library) or `libpetibm.so` (shared library) to link against.
The library is located in `$PETIBM_INSTALL/lib`; use `-L$PETIBM_INSTALL/lib -lpetibm` to link.

There may be other libraries in `$PETIBM_INSTALL/lib`, depending on how the users configure PetIBM.
For example, if the users use the configuration flag `--enable-yamlcpp`, then, the yaml-cpp library will be installed along with the PetIBM library in `$PETIBM_INSTALL/lib`.
(Note: on some systems or OS, third-party libraries may be installed in `$PETIBM_INSTALL/lib64`.)

Also, linking against PETSc library is necessary.
Flags for linking against PETSc can be found in the file `$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables`.
(Look for the line starting with `PETSC_WITH_EXTERNAL_LIB = `.)

Alternatively, users can let the PetIBM library automatically handle dependencies with the `libtool` program.
PetIBM is build using `libtool`, which records dependencies information in an abstract layer `$PETIBM_INSTALL/lib/libpetibm.la`.
Although users can directly look into the `.la` file, the proper way to use `libtool` is:

```
libtool --tag=CXX --mode=link $CXX -o <executable name> <object files> \
    $LINKER_FLAGS $PETIBM_INSTALL/lib/libpetibm.la
```

or ,if compiling and linking at once,

```
libtool --tag=CXX --mode=link $CXX -o <executable name> <source files>  \
    $CXX_FLAGS $LINKER_FLAGS $PATH_TO_HEADERS $PETIBM_INSTALL/lib/libpetibm.la
```

