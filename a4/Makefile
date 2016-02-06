CONFIG ?= config
include $(CONFIG).mk

BDIR = ./bin
SDIR = ./src
IDIR = ./inc
ODIR = ./obj


#===============================
#	  PETSC NUMERICAL LIB 
#===============================

PETSC_INC = -I $(PETSC_DIR)/include
PETSC_LIB = -L $(PETSC_DIR)/lib

#===============================
#	  XERCES C++ XML PARSING 
#===============================

# XERCES_INC = -I $(XERCES_DIR)/include
# XERCES_LIBS += -L $(XERCES_DIR)/lib

#===============================
#			PUMI 
#===============================

PUMI_INC = -I $(PUMI_DIR)/include
PUMI_LIBS += -L $(PUMI_DIR)/lib

#===============================
#		GOOGLE TEST 
#===============================

# Where to find user code.
USER_DIR = $(GTEST_DIR)/samples

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -isystem $(GTEST_DIR)

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread --pedantic-errors

#now that we are using mpich, clang errors on long long
MPIFLAGS += -Wall -g -I $(IDIR) --pedantic-errors

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = $(BDIR)/app_test $(BDIR)/convergence_test

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
				$(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(TESTS)

#there are other things in lib, so only remove libraries explictly by name
clean :
	rm -vf $(TESTS) $(BDIR)/* $(LDIR)/gtest.a $(LDIR)/gtest_main.a $(ODIR)/*.o

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(ODIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest-all.cc -o $@

$(ODIR)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest_main.cc -o $@

$(LDIR)/gtest.a : $(ODIR)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(LDIR)/gtest_main.a : $(ODIR)/gtest-all.o $(ODIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.


$(ODIR)/MeshAdjReorder.o : 	$(SDIR)/MeshAdjReorder.cc \
					  		$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/MeshBuilder_unittest.o : 	$(SDIR)/MeshBuilder_unittest.cc \
					  				$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/MeshBuilder.o : $(SDIR)/MeshBuilder.cc $(IDIR)/MeshBuilder.h \
					  	$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/SimpleRectMesh_unittest.o : $(SDIR)/SimpleRectMesh_unittest.cc \
					  				$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/ElasticAnalysis2D.o : 	$(SDIR)/ElasticAnalysis2D.cc \
								$(IDIR)/ElasticAnalysis2D.h \
								$(IDIR)/RecoverAtIntegrationPoints.h
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/ElasticAnalysis_unittest.o : 	$(SDIR)/ElasticAnalysis_unittest.cc \
										$(IDIR)/ElasticAnalysis2D.h \
										$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/ShapeFunctionOrdering_unittest.o : 	$(SDIR)/ShapeFunctionOrdering_unittest.cc \
										$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@


$(ODIR)/StiffnessContributor2D.o : 	$(SDIR)/StiffnessContributor2D.cc \
								$(IDIR)/StiffnessContributor2D.h \
								$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/ForceContributor2D.o : 	$(SDIR)/ForceContributor2D.cc \
								$(IDIR)/ForceContributor2D.h
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/AlgebraicSystem.o : $(SDIR)/AlgebraicSystem.cc \
							$(IDIR)/AlgebraicSystem.h
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/AlgebraicSystem_unittest.o : 	$(SDIR)/AlgebraicSystem_unittest.cc \
										$(IDIR)/AlgebraicSystem.h \
										$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/GeometryMappings.o : 	$(SDIR)/GeometryMappings.cc \
								$(IDIR)/GeometryMappings.h \
								$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/GeometryMappings_unittest.o : 	$(SDIR)/GeometryMappings_unittest.cc \
										$(IDIR)/GeometryMappings.h \
										$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/StiffnessContributor2D_unittest.o : $(SDIR)/StiffnessContributor2D_unittest.cc \
							$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(ODIR)/Convergence_test.o : 	$(SDIR)/Convergence_test.cc \
								$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) $(PETSC_INC) -c $< -o $@

#pumi library is no longer rebuilt every time
$(ODIR)/mpi_and_petsc_gtest_main.o : $(SDIR)/mpi_and_petsc_gtest_main.cc $(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@


$(BDIR)/app_test : 	$(ODIR)/mpi_and_petsc_gtest_main.o $(ODIR)/MeshBuilder.o $(ODIR)/MeshBuilder_unittest.o \
				$(ODIR)/SimpleRectMesh_unittest.o \
				$(ODIR)/MeshAdjReorder.o \
				$(ODIR)/ElasticAnalysis_unittest.o \
				$(ODIR)/ElasticAnalysis2D.o \
				$(ODIR)/ShapeFunctionOrdering_unittest.o \
				$(ODIR)/StiffnessContributor2D.o \
				$(ODIR)/ForceContributor2D.o \
				$(ODIR)/AlgebraicSystem.o \
				$(ODIR)/AlgebraicSystem_unittest.o \
				$(ODIR)/GeometryMappings.o \
				$(ODIR)/GeometryMappings_unittest.o \
				$(ODIR)/StiffnessContributor2D_unittest.o \
				$(LDIR)/gtest.a
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(PETSC_LIB) $(PUMI_LIBS) $^ -lpetsc -lmds -lapf -lgmi -lpcu -llion -lmth -lpthread  -o $@


$(BDIR)/convergence_test : 	$(ODIR)/mpi_and_petsc_gtest_main.o \
							$(ODIR)/Convergence_test.o \
							$(ODIR)/MeshBuilder.o \
							$(ODIR)/MeshAdjReorder.o \
							$(ODIR)/ElasticAnalysis2D.o \
							$(ODIR)/StiffnessContributor2D.o \
							$(ODIR)/ForceContributor2D.o \
							$(ODIR)/AlgebraicSystem.o \
							$(ODIR)/GeometryMappings.o \
							$(LDIR)/gtest.a
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(PETSC_LIB) $(PUMI_LIBS) $^ -lpetsc -lmds -lapf -lgmi -lpcu -llion -lmth -lpthread -o $@

#===============================================

german: $(BDIR)/bug.out

$(ODIR)/bug.o : bug.cc
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PETSC_INC) $(PUMI_INC) -c $< -o $@

$(BDIR)/bug.out : $(ODIR)/bug.o $(ODIR)/MeshBuilder.o
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(PETSC_LIB) $(PUMI_LIBS) $^ -lpetsc -lmds -lapf -lgmi -lpcu -llion -lmth -lpthread  -o $@


Goofy : $(BDIR)/test_xml_reader.out

$(ODIR)/test_xml_reader.o : $(SDIR)/test_xml_reader.cc
	mpicxx $(CPPFLAGS) $(CXXFLAGS) -I $(IDIR) $(PUMI_INC) $(XERCES_INC) -c $< -o $@

$(ODIR)/vtkSAX2Handler.o : $(SDIR)/vtkSAX2Handler.cc
	mpicxx $(CPPFLAGS) $(CXXFLAGS) -I $(IDIR) $(PUMI_INC) $(XERCES_INC) -c $< -o $@

$(BDIR)/test_xml_reader.out : $(ODIR)/test_xml_reader.o $(ODIR)/vtkSAX2Handler.o
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $(XERCES_LIBS) $^\
			-lpthread -lxerces-c -o $@

#old class not used anymore
$(ODIR)/BandedSymmetricMatrix.o : 	$(SDIR)/BandedSymmetricMatrix.cc \
								$(IDIR)/BandedSymmetricMatrix.h \
								$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@

$(ODIR)/BandedSymmetricMatrix_unittest.o : 	$(SDIR)/BandedSymmetricMatrix_unittest.cc \
								$(GTEST_HEADERS)
	mpicxx $(CPPFLAGS) $(MPIFLAGS) $(PUMI_INC) -c $< -o $@