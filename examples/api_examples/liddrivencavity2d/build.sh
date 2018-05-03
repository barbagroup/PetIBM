#!/bin/sh


# ============================================================================
# Edit the following variables
# ============================================================================

# The path to PETSc installation
PETSC_DIR=

# The arch of specific PETSc
PETSC_ARCH=

# The path to PetIBM installation
PETIBM_DIR=

# Extra flags. Please provide the include path and linker flags to the libraries
# that were not build through PetIBM configuration. For example, if you built
# and install yaml-cpp by yourself.
EXTRA_INCLUDES=


# ============================================================================
# The following variables are automatically obtained
# ============================================================================

# path to the file of PETSc varaibles
PETSCVARIABLES=${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables

# the compiler used to build PETSc
CXX=" \
    `grep "CXX =" $PETSCVARIABLES | \
    sed -e 's/CXX\ =\ //' -e 's/^[ \t]*//'`"

# manually set compilation flags for debuging purpose
CXX_FLAGS="--std=c++11 -g"

# the include path for PETSc
PETSC_CC_INCLUDES=" \
    `grep "PETSC_CC_INCLUDES =" $PETSCVARIABLES \
    | sed -e 's/.*=//' -e 's/^[ \t]*//'`"

# the include path for PetIBM
PETIBM_INCLUDES=-I$PETIBM_DIR/include


# ============================================================================
# Output information
# ============================================================================

echo "[Compiler]"
echo $CXX
echo
echo "[C++ FLAGS]"
echo $CXX_FLAGS
echo
echo "[PETSc INCLUDE FLAGS]"
echo $PETSC_CC_INCLUDES
echo
echo "[PetIBM INCLUDE FLAGS]"
echo $PETIBM_INCLUDES
echo

# ============================================================================
# Start to build
# ============================================================================

echo "Start building:" 
libtool --tag=CXX --mode=link $CXX -o a.out main.cpp \
    $CXX_FLAGS $PETSC_CC_INCLUDES \
    $PETIBM_INCLUDES $EXTRA_INCLUDES \
    $PETIBM_DIR/lib/libpetibm.la
