CC=gcc

ifeq ($(CC), gcc)
	MPIFLAGS+=--std=c++0x -Wno-long-long
else
	MPIFLAGS+=-Wno-c++11-long-long
endif

# PETSC_DIR=/tmp/packages/petsc-3.6.1
# LDIR=$(MAIN_DIR)/lib
# PUMI_DIR=$(LDIR)/pumi
# GTEST_DIR=/tmp/googletest


# MPIFLAGS+="-cxx=clang++ -Wno-c++11-long-long"

PETSC_DIR = /home/shivaebola/Documents/Software/petsc/petsc-3.6.1

LDIR = /home/shivaebola/GitLab/FEP/lib
XERCES_DIR = $(LDIR)/xerces-c-3.1.2/build
PUMI_DIR = $(LDIR)/pumi
# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR =  ~/GitLab/googletest

