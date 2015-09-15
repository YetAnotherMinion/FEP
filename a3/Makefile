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

#use DEPS to require a rebuild whenever a header file changes
_DEPS = 
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than main target
_OBJ = 
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread --pedantic-errors
#now that we are using mpich, clang errors on long long
MPIFLAGS += -Wall -g -I $(IDIR) --pedantic-errors

all: $(BDIR)/a3

$(BDIR)/a3 : $(ODIR)/a3.o $(OBJ)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(CXXFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -lpthread -o $@ 

$(ODIR)/a3.o : $(SDIR)/a3.cc 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(CXXFLAGS) $(PUMI_INC) -c $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/* $(ODIR)/*.o