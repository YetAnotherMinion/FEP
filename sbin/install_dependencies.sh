#!/bin/bash

set -eux
self=$(readlink -f $(dirname "$0"))
export PATH=${self}:$PATH


# Default to only download and install missing dependencies
# This way devs can modify the git folders
clean_build=false

# prefix is where we will install everything
prefix=$HOME/tmp

while getopts ":p:c" opt; do
    case $opt in
        c)
            clean_build=true
            ;;
        p)
            prefix=$OPTARG
            ;;
        :)
            echo "Option $OPTARG requires an argument" >&2
            exit 1
    esac
done

#### Configure Your Install Here ####

# override compiler defaults here
CC=${CC:-gcc}
CXX=${CXX:-g++}

# By default I use all but one core for compiling
max_make_threads=$(expr $(getconf _NPROCESSORS_ONLN) - 1)
max_make_load=${max_make_threads}

set +u
export PATH=$PATH:${prefix}/bin
export CPATH=$CPATH:${prefix}/inc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${prefix}/lib
set -u
package_dir=${prefix}/src

# Download locations for dependencies can be changed here to use internal
# artifact caches
MPICH_DL_LOC="${MPICH_DL_LOC:-http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz}"
BLAS_LAPACK_DL_LOC="${BLAS_LAPACK_DL_LOC:-http://ftp.mcs.anl.gov/pub/petsc/externalpackages/fblaslapack-3.4.2.tar.gz}"
PUMI_TEST_MESHES_DL_LOC="${PUMI_TEST_MESHES_DL_LOC:-https://www.scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz}"

mpich_tarball=${package_dir}/mpich-3.2.tar.gz
mpich_chksum="083c51655b4355827bd7fa4fe528046e2bc77b7747d869ff87b79fa324c3cc2a9b5640ccb7271490ccc0dd627e354a33a449bbab448501bbfddcfe5f999ee717"

blas_tarball=${package_dir}/fblaslapack-3.4.2.tar.gz
blas_chksum="81020ed02a2f7bd467f8417ff433c5a809a5fa49bbc77ae34e2e157b36aecaa49ea936faf6599465b0df7e97d669a267087d26e8cdd94ba7d80395c1c3380897"

pumi_test_meshes_tarball=${package_dir}/pumi_test_meshes.tar.gz
pumi_test_meshes_chksum="5a1530b27eb4cc88958f5758a64f10eadd0f564dc0959bd3e10a2d51721eb41dff529734cbe28d9e83ea5068d47c56952979c3c845870f18aecb95900f0cb5b8"

gtest_dir=${package_dir}/googletest
petsc_dir=${package_dir}/petsc
superlu_dir=${package_dir}/git.superlu
sowing_dir=${package_dir}/git.sowing
core_dir=${package_dir}/core
test_meshes_dir=${package_dir}/test_meshes

#### End of user configuration ####


if $clean_build; then
    rm -rf $mpich_tarball $blas_tarball $pumi_test_meshes_tarball
    rm -rf $gtest_dir $petsc_dir $core_dir $test_meshes_dir
fi

mkdir -p $package_dir

# Do all network traffic at one point in the script in order to increase
# the chances of a thundering herd
memoized_curl ${MPICH_DL_LOC} $mpich_tarball $mpich_chksum
memoized_curl ${BLAS_LAPACK_DL_LOC} $blas_tarball $blas_chksum
memoized_curl ${PUMI_TEST_MESHES_DL_LOC} $pumi_test_meshes_tarball $pumi_test_meshes_chksum


if [[ ! -d $petsc_dir ]]; then
    git clone -b maint https://bitbucket.org/petsc/petsc $petsc_dir
fi

if [[ ! -d $gtest_dir ]]; then
    git clone https://github.com/YetAnotherMinion/googletest.git $gtest_dir
fi

if [[ ! -d $core_dir ]]; then
    git clone https://github.com/SCOREC/core.git $core_dir
fi

if [[ ! -d $test_meshes_dir ]]; then
    git clone https://github.com/SCOREC/fep.git $test_meshes_dir
fi

if [[ ! -d $superlu_dir ]]; then
    git clone https://github.com/xiaoyeli/superlu $superlu_dir
fi

if [[ ! -d $sowing_dir ]]; then
    git clone https://bitbucket.org/petsc/pkg-sowing.git $sowing_dir
fi

cd ${package_dir}
#### Build MPICH2 ####

tar -xzf $mpich_tarball
cd mpich-3.2
./configure --prefix=${prefix}
make -j ${max_make_threads} -l ${max_make_load}
make install

#### Build PETSC ####
cd $petsc_dir
export PETSC_ARCH=linux-mpich
export PETSC_DIR=${petsc_dir}

./configure --prefix=$prefix \
    --download-fblaslapack=${blas_tarball} \
    --download-superlu=${superlu_dir} \
    --download-sowing=${sowing_dir} \
    --with-pthread=yes \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --CXX=$CXX \
    --CC=$CC

make -j ${max_make_threads} -l ${max_make_load}  all
make install
# test the installed library
make test

#### Build PUMI ####
mkdir -p ${core_dir}/build
cd ${core_dir}/build
tar -xzf ${pumi_test_meshes_tarball} -C ${core_dir}

# right now we only support using MPICH
export MPICH_CXX=$CXX
export MPICH_CC=$CC

# need to copy the script into a specific location because it calls
# subshells that use relative paths. We are hotpatching the regular
# install process anyway
cp $(which custom_pumi_configure.sh) ${core_dir}/custom_pumi_configure.sh
source ../custom_pumi_configure.sh ${prefix}/lib/pumi ${package_dir}/core/meshes
make -j ${max_make_threads} -l ${max_make_load}
make install -k
ctest -W

