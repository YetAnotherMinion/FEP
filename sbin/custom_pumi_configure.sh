#!/bin/bash

#first arg is installation directory
#second arg is test mesh location directory
if [[ ! -z ${MPI_VENDOR-} ]]; then
    case $MPI in
	    "openmpi")
            export OMPI_CXX=$CXX
            export OMPI_CC=$CC
            ;;
	    "mpich2")
            export MPICH_CXX=$CXX
            export MPICH_CC=$CC
            ;;
    esac
fi

cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
  -DENABLE_THREADS=OFF \
  -DENABLE_ZOLTAN=OFF \
  -DCMAKE_INSTALL_PREFIX="$1" \
  -DIS_TESTING=True \
  -DMESHES="$2"
