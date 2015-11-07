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

CFLAGS = -Wall -g -I $(IDIR)  $(PUMI_INC) --pedantic-errors -lmds -lapf -lgmi -lpcu  

#use DEPS to require a rebuild whenever a header file changes
_DEPS = 
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

#hold all the objects other than server
#this is so that whenever something changes the server is recompiled
_OBJ = 
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

#now that we are using mpich, clang errors on long long
MPIFLAGS += -Wall -g -I $(IDIR) --pedantic-errors

all: $(BDIR)/a2

$(BDIR)/a2 : $(ODIR)/a2.o $(OBJ)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_LIBS) $< -lmds -lapf -lgmi -lpcu -llion -lmth -lpthread -o $@

$(ODIR)/a2.o : $(SDIR)/a2.cc 
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@ 

#prevent anything from happening to a file named clean
.FAKE : clean

clean:
	$(RM) $(BDIR)/* $(ODIR)/*.o
