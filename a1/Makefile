CONFIG ?= config
include $(CONFIG).mk

BDIR = ./bin
SDIR = ./src
IDIR = ./inc
ODIR = ./obj

#===============================
#			PUMI 
#===============================

PUMI_INC = -I $(PUMI_DIR)/include
PUMI_LIBS += -L $(PUMI_DIR)/lib

#--verbose adds a lot of information to output and expands on include search path
#-L$(Lgmi)/libgmi  -L$(Lapf)/libapf  -L$(Lpcu)/libpcu  -L$(Lmds)/libmds


CFLAGS = -Wall -g -I $(IDIR) --pedantic-errors
#LIBS =  $(Lapf)/libapf.a $(Lgmi)/libgmi.a $(Lmds)/libmds.a $(Lpcu)/libpcu.a

#use DEPS to require a rebuild whenever a header file changes
_DEPS = 
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than server
#this is so that whenever something changes the server is recompiled
_OBJ = 
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

#`pkg-config --cflags --libs libmds`

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread --pedantic-errors

#now that we are using mpich, clang errors on long long
MPIFLAGS += -Wall -g -I $(IDIR) --pedantic-errors

all: ./bin/a1tet ./bin/a1quad ./bin/a1mix

$(BDIR)/a1tet : $(ODIR)/a1tet.o
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(LIBS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1tet.o : $(SDIR)/a1tet.cc 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@ 

$(BDIR)/a1quad : $(ODIR)/a1quad.o 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1quad.o : $(SDIR)/a1quad.cc 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@ 

$(BDIR)/a1mix : $(ODIR)/a1mix.o
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a1mix.o : $(SDIR)/a1mix.cc 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/* $(ODIR)/*.o
