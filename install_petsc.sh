#! /usr/bin/env sh

#takes three arguments
# first is install directory
# second is the full path to petsc directory
# third is full path to directory where packages have been downloaded without
# trailing slash

./configure --prefix=$1 \
  PETSC_ARCH=linux-mpich  PETSC_DIR=$2 \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --CXX=clang++ --CC=clang \
  --download-fblaslapack=$3/fblaslapack-3.4.2.tar.gz \
  --download-superlu=$3/superlu_4.3.tar.gz \
  --with-pthread=yes
