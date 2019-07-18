#!/bin/bash

# Install PETSc-3.11.3 in debug mode.
# Use the cached PETSc if mpicxx is detected.

PETSC_DIR=$1
PETSC_ARCH=$2

if [ -f "$PETSC_DIR/$PETSC_ARCH/bin/mpicxx" ]; then
  echo "Using cached PETSc"
  echo "Configuring PETSc"
  cd $PETSC_DIR
  ./configure PETSC_ARCH=$PETSC_ARCH \
    --with-cc=$CC --with-cxx=$CXX --with-fc=0 \
    --COPTFLAGS="-O0" --CXXOPTFLAGS="-O0" --FOPTFLAGS="-O0" \
    --with-debugging=1 \
    --download-f2cblaslapack --download-hdf5 --download-openmpi
else
  echo "Downloading PETSc"
  URL="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.3.tar.gz"
  TARBALL=/tmp/petsc-lite-3.11.3.tar.gz
  wget $URL -O $TARBALL
  tar xfz $TARBALL -C $PETSC_DIR --strip-components=1
  rm -f $TARBALL
  echo "Configuring and building PETSc"
  cd $PETSC_DIR
  ./configure PETSC_ARCH=$PETSC_ARCH \
    --with-cc=$CC --with-cxx=$CXX --with-fc=0 \
    --COPTFLAGS="-O0" --CXXOPTFLAGS="-O0" --FOPTFLAGS="-O0" \
    --with-debugging=1 \
    --download-f2cblaslapack --download-hdf5 --download-openmpi
  make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
  make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH test
fi
