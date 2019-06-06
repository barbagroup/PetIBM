# Dockerfile for PetIBM-0.4.1

FROM nvidia/cuda:8.0-devel-ubuntu16.04
MAINTAINER Olivier Mesnard <mesnardo@gwu.edu>

# Install base system.
COPY ssh_config /root/.ssh/config
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        autotools-dev \
        gfortran \
        cmake \
        git \
        pkg-config \
        wget \
        curl \
        unzip && \
    # OpenMPI
    apt-get install -y --no-install-recommends \
        libopenmpi-dev \
        libopenmpi1.10 \
        openmpi-bin \
        openmpi-common \
        openmpi-doc && \
    # For PETSc
    apt-get install -y --no-install-recommends \
        flex \
        bison \
        python-dev && \
    # For RDMA/InfiniBand
    apt-get install -y --no-install-recommends \
        libibverbs1 ibverbs-utils librdmacm1 rdmacm-utils \
        ibutils ibacm libibcm1 libibmad5 libibumad3 opensm libopensm5a \
        srptools perftest infiniband-diags ibsim-utils \
        libmthca1 libmlx4-1 libmlx5-1 mstflint \
        libdapl2 dapl2-utils && \
    # configure ssh server and keys
    apt-get install -y --no-install-recommends \
        openssh-server \
        openssh-client && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /var/run/sshd && \
    ssh-keygen -A && \
    sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config && \
    sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd && \
    ssh-keygen -f /root/.ssh/id_rsa -t rsa -N '' && \
    chmod 600 /root/.ssh/config && \
    chmod 700 /root/.ssh && \
    cp /root/.ssh/id_rsa.pub /root/.ssh/authorized_keys

# Build PETSc-3.10.2 in optimized mode.
RUN VERSION=3.10.2 && \
    TARBALL=petsc-lite-${VERSION}.tar.gz && \
    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${TARBALL} -P /tmp && \
    PETSC_DIR=/opt/petsc/${VERSION} && \
    PETSC_ARCH=linux-gnu-openmpi-opt && \
    mkdir -p ${PETSC_DIR} && \
    tar -xzf /tmp/${TARBALL} -C ${PETSC_DIR} --strip-components=1 && \
    cd ${PETSC_DIR} && \
    ./configure --PETSC_DIR=${PETSC_DIR} --PETSC_ARCH=${PETSC_ARCH} \
      --with-cc=mpicc \
      --with-cxx=mpicxx \
      --with-fc=mpif90 \
      --COPTFLAGS=-O3 \
      --CXXFLAGS=-O3 \
      --FOPTFLAGS=-O3 \
      --with-debugging=0 \
      --download-hdf5 \
      --download-fblaslapack \
      --download-hypre \
      --download-ptscotch \
      --download-metis \
      --download-parmetis \
      --download-superlu_dist && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all && \
    rm -f /tmp/${TARBALL}

# Build AmgX-2.0.
RUN VERSION=2.0 && \
    TARBALL=master.tar.gz && \
    wget https://github.com/NVIDIA/AMGX/archive/${TARBALL} -P /tmp && \
    SRCDIR=/opt/amgx/${VERSION} && \
    BUILDDIR=${SRCDIR}/build && \
    mkdir -p ${SRCDIR} ${BUILDDIR} && \
    tar -xzf /tmp/${TARBALL} -C ${SRCDIR} --strip-components=1 && \
    cd ${BUILDDIR} && \
    cmake ${SRCDIR} \
      -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_INSTALL_PREFIX="/usr/local" \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS_PROFILE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_CXX_FLAGS_PROFILE="-O3 -DNDEBUG" \
      -DMPI_CXX_COMPILER=mpicxx \
      -DMPI_CXX_COMPILE_FLAGS="-O3" \
      -DMPI_C_COMPILER=mpicc \
      -DMPI_C_COMPILE_FLAGS="-O3" \
      -DCUDA_ARCH="35 37" \
      -DCUDA_HOST_COMPILER=/usr/bin/mpicc && \
    make -j"$(nproc)" all && \
    make install && \
    rm -f /tmp/${TARBALL} /opt/amgx/srcTarball.txt && \
    ldconfig ${SRCDIR}/lib

ENV PETSC_DIR=/opt/petsc/3.10.2
ENV PETSC_ARCH=linux-gnu-openmpi-opt

# Build and install PetIBM-0.4.1.
RUN VERSION=0.4.1 && \
    TARBALL=v${VERSION}.tar.gz && \
    wget https://github.com/barbagroup/PetIBM/archive/${TARBALL} -P /tmp && \
    SRCDIR=/opt/petibm/${VERSION} && \
    BUILDDIR=${SRCDIR}/build && \
    mkdir -p ${SRCDIR} ${BUILDDIR} && \
    tar -xzf /tmp/${TARBALL} -C ${SRCDIR} --strip-components=1 && \
    cd ${BUILDDIR} && \
    ${SRCDIR}/configure --prefix=/usr/local \
      CXX=mpicxx \
      CXXFLAGS="-O3 -w -std=c++14" \
      --enable-static=no \
      --with-petsc-dir=${PETSC_DIR} \
      --with-petsc-arch=${PETSC_ARCH} \
      --with-cuda-dir=/usr/local/cuda-8.0 \
      --with-amgx-dir=/usr/local \
      --enable-amgxwrapper \
      --enable-yamlcpp \
      --enable-gtest && \
    make -j"$(nproc)" all && \
    make check && \
    make install && \
    ldconfig /usr/local/lib && \
    rm -f /tmp/${TARBALL}

EXPOSE 23
CMD ["/usr/sbin/sshd", "-D", "-p", "23"]
