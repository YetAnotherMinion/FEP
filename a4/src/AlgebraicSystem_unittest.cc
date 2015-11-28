#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <apf.h>
#include <apfDynamicMatrix.h>

#include <gtest/gtest.h>

#include <petscmat.h>
#include <petscksp.h>

#include "AlgebraicSystem.h"

class AlgebraicSystemTest : public testing::Test
{
protected:
	std::size_t nGlobalDOFs;


	apf::DynamicMatrix ke_1;
	apf::NewArray<int> mapping_1;
	uint32_t ke_1_size;
	/*matrix designed to fit on lower diagonal of global matrix*/
	apf::DynamicMatrix ke_2;
	apf::NewArray<int> mapping_2;
	uint32_t ke_2_size;
	/*fixture for some displacements testing*/
	std::vector<double> local_disp;
	std::vector<uint64_t> local_mapping;


	virtual void SetUp() {
		nGlobalDOFs = 12; /*must pick avalue above 4*/
		/*================================*/
		ke_1_size = 4;
		mapping_1.allocate(ke_1_size);	
		ke_1.setSize(ke_1_size, ke_1_size);
		for(std::size_t ii = 0; ii < ke_1_size; ++ii) {
			for(std::size_t jj = 0; jj < ke_1_size; ++jj) {
				if(ii == jj) {
					ke_1(ii, jj) = ii * jj;
				} else {
					ke_1(ii, jj) = 0.0;
				}
			}
		}
		/*this structure will spread out the submatrix
		* all over the larger matrix*/
		mapping_1[0] = 1;
		mapping_1[1] = 3;
		mapping_1[2] = nGlobalDOFs-3;
		mapping_1[3] = nGlobalDOFs-1;
		/*==================================*/
		ke_2_size = 4;
		mapping_2.allocate(ke_2_size);	
		ke_2.setSize(ke_2_size, ke_2_size);
		for(std::size_t ii = 0; ii < ke_2_size; ++ii) {
			for(std::size_t jj = 0; jj < ke_2_size; ++jj) {
				if(ii == jj) {
					/*on the main diagonal*/
					ke_2(ii, jj) = 10.0 * ii + jj + 1;
				} else if(ii > jj) {
					ke_2(ii,jj) = 7.0;
				} else {
					ke_2(ii, jj) = 0.0;
				}
			}
		}
		/*this structure will spread out the submatrix
		* all over the larger matrix*/
		mapping_2[0] = nGlobalDOFs-4;
		mapping_2[1] = nGlobalDOFs-3;
		mapping_2[2] = nGlobalDOFs-2;
		mapping_2[3] = nGlobalDOFs-1;
		/*================================*/
		local_disp.clear();
		local_disp.push_back(6.0);
		local_disp.push_back(8.0);
		local_disp.push_back(7.0);
		local_mapping.clear();
		local_mapping.push_back(nGlobalDOFs-4);
		local_mapping.push_back(nGlobalDOFs-1);
		local_mapping.push_back(0);
	}

	virtual void TearDown() {

	}

};

TEST_F(AlgebraicSystemTest, AssembleInsertionTest) {
	/*in this test we check that insertions of stiffness matrix
	* are additive, that is if we insert the same local stiffness 
	* matrix into the global matrix that the contributions are
	* additive. */

	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(nGlobalDOFs);

	linsys.beginAssembly();

	uint32_t nReps = 4;
	for(std::size_t ii = 0; ii < nReps; ++ii) {
		linsys.assemble(ke_1, mapping_1, ke_1_size);
	}
	linsys.synchronize();

	PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;

	PetscInt row;
	/*check that the rows above the diagonal are correct*/
	for(std::size_t ii = 0; ii < ke_1_size; ++ii){
		row = mapping_1[ii];
		/*we only expect to see the upper triangular portion*/
		MatGetRowUpperTriangular(linsys.K);
		MatGetRow(linsys.K, row, &ncols, &cols, &vals);
		/*cross check with the initial submatrix*/
		assert(ncols >= 0);
		for(std::size_t jj = 0; jj < static_cast<std::size_t>(ncols); ++jj) {
			std::size_t kk = -1;
			/*find the mapping for this column back to local
			* this is done by looping over every mapping until we find it
			* this is hacky but ok since the size is only 4*/
			for(std::size_t ll = 0; ll < ke_1_size; ++ll) {
				if(mapping_1[ll] == cols[jj]) {
					kk = ll;
					break; /*stop once we have found a mapping*/
				}
			}
			if(kk != static_cast<std::size_t>(-1)) {
				EXPECT_TRUE(kk < ke_1_size);
				EXPECT_TRUE(ii < ke_1_size);
				/*the expected value is nReps times the local matrix*/
				EXPECT_FLOAT_EQ((ke_1(ii, kk) * nReps), vals[jj]);
				
			} else {
				/*we should always find mapping*/
				FAIL();
			}
		}
		MatRestoreRow(linsys.K, row, &ncols, &cols, &vals);
		MatRestoreRowUpperTriangular(linsys.K);
	}

}

TEST_F(AlgebraicSystemTest, IgnoreLowerDiagonals) {
	/*In this test we insert a local stiffness matrix
	* that has components above and below the main diagonal
	* of the global matrix, we ensure that these lower
	* contributions are not added to the upper triangular
	* portion of the global matrix by the underlying 
	* implementation. Basically we are checking that the
	* library provided matrix behaves the way we think
	* we are using it. */


	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);

	linsys.beginAssembly();
	linsys.assemble(ke_2, mapping_2, ke_2_size);
	linsys.synchronize();

	PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;

	PetscInt row;
	/*check that the rows above the diagonal are correct*/
	for(std::size_t ii = 0; ii < ke_2_size; ++ii){
		row = mapping_2[ii];
		/*we only expect to see the upper triangular portion*/
		MatGetRowUpperTriangular(linsys.K);
		MatGetRow(linsys.K, row, &ncols, &cols, &vals);
		/*cross check with the initial submatrix*/
		assert(ncols >= 0);
		for(std::size_t jj = 0; jj < static_cast<std::size_t>(ncols); ++jj) {
			std::size_t kk = -1;
			/*find the mapping for this column back to local
			* this is done by looping over every mapping until we find it
			* this is hacky but ok since the size is only 4*/
			for(std::size_t ll = 0; ll < ke_2_size; ++ll) {
				if(mapping_2[ll] == cols[jj]) {
					kk = ll;
					break; /*stop once we have found a mapping*/
				}
			}
			if(kk != static_cast<std::size_t>(-1)) {
				EXPECT_TRUE(kk < ke_2_size);
				EXPECT_TRUE(ii < ke_2_size);
				/*only test when */
				EXPECT_FLOAT_EQ(ke_2(ii, kk), vals[jj]);
				
			} else {
				/*we should always find mapping*/
				FAIL();
			}
		}
		MatRestoreRow(linsys.K, row, &ncols, &cols, &vals);
		MatRestoreRowUpperTriangular(linsys.K);
	}
}

TEST_F(AlgebraicSystemTest, PreventPrematureAssembly) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*try to call assemble before calling beginAssembly
	* here we test for failure in form of execption*/
	EXPECT_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size), std::invalid_argument);
}

TEST_F(AlgebraicSystemTest, AllowCorrectAssembly) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*try to call assemble after calling beginAssembly
	* here we do not want an exception*/
	linsys.beginAssembly();
	EXPECT_NO_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size));
}

TEST_F(AlgebraicSystemTest, PreventMutltipleAssembly) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	EXPECT_NO_THROW(linsys.beginAssembly());
	EXPECT_THROW(linsys.beginAssembly(), std::logic_error);
}

TEST_F(AlgebraicSystemTest, RepeatedConstraints) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
	/*over constraining the system is allowed, even if we add the same entity
	* twice, it is expected to overwrite the exisiting displacements. So to 
	* verify we change one of the displacements*/
	this->local_disp[0] += 27.0;
	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
	/*the known displacements are private, so we must go ahead and solve
	* the empty system to then extract the solution vector in order to check
	* if the overwriting was sucessful. We are only solving a tiny system,
	* so this behaviour is acceptable*/
	linsys.beginAssembly();
	linsys.synchronize();
	linsys.solve();

	std::vector<double> solution;
	solution.clear();
	linsys.extractDisplacement(solution);
	EXPECT_EQ(nGlobalDOFs, solution.size());
	for(std::size_t ii = 0; ii < nGlobalDOFs; ++ii) {
		if(ii == this->local_mapping[0]) {
			EXPECT_EQ(this->local_disp[0], solution[ii]);
		} else if(ii == this->local_mapping[1]) {
			EXPECT_EQ(this->local_disp[1], solution[ii]);
		} else if(ii == this->local_mapping[2]) {
			EXPECT_EQ(this->local_disp[2], solution[ii]);
		} else {
			EXPECT_FLOAT_EQ(0.0, solution[ii]);
		}
	}
}

TEST_F(AlgebraicSystemTest, OutOfRangeContraint) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*corrupt by making this one clearly out of range*/
	this->local_mapping[1] = nGlobalDOFs+1;

	EXPECT_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping), std::invalid_argument);
}

TEST_F(AlgebraicSystemTest, ExtraMappingsAreIgnored) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*this one is clearly out of range but it will not be used since there are 
	* fewer displacements than there are mappings, so this last one will not be used*/
	local_mapping.push_back(nGlobalDOFs+1);
	/*quick test for consitency, this should always pass, only included
	* to prevent others from dicking with the test fixture and messing
	* this up*/
	ASSERT_EQ(this->local_mapping.size(), this->local_disp.size()+1);

	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
}

TEST_F(AlgebraicSystemTest, ExtractSolutionWithConstraints) {
	AlgebraicSystem linsys(this->nGlobalDOFs);

	linsys.addBoundaryConstraint(this->local_disp, this->local_mapping);
	linsys.beginAssembly();
	linsys.synchronize();
	linsys.solve();

	std::vector<double> solution;
	solution.clear();
	linsys.extractDisplacement(solution);
	EXPECT_EQ(nGlobalDOFs, solution.size());
	for(std::size_t ii = 0; ii < nGlobalDOFs; ++ii) {
		if(ii == this->local_mapping[0]) {
			EXPECT_EQ(this->local_disp[0], solution[ii]);
		} else if(ii == this->local_mapping[1]) {
			EXPECT_EQ(this->local_disp[1], solution[ii]);
		} else if(ii == this->local_mapping[2]) {
			EXPECT_EQ(this->local_disp[2], solution[ii]);
		} else {
			EXPECT_FLOAT_EQ(0.0, solution[ii]);
		}
	}
}

	// PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	// MatView(linsys.K, PETSC_VIEWER_STDOUT_WORLD);
	// VecView(linsys.F, PETSC_VIEWER_STDOUT_WORLD);

#define MATRIX_SIZE 20

class MatrixMathTest : public testing::Test
{
protected:
	PetscScalar A[MATRIX_SIZE*MATRIX_SIZE];
	PetscScalar b[MATRIX_SIZE];
	PetscScalar x[MATRIX_SIZE];
	PetscInt rows[MATRIX_SIZE];
	Mat AA;
	Vec bb;
	Vec xx;
	KSP solver;
	apf::DynamicMatrix ke;
	apf::NewArray<int> ke_mapping;

	virtual void SetUp() {
		/*obtained from matlab by:  round(gallery('randcorr',20).* 10^5))./10^5 to reduce the text length of the floats*/
		PetscScalar tmp_A[MATRIX_SIZE*MATRIX_SIZE] = {1,0.11042,0.15975,0.07215,-0.16194,0.12573,0.16026,0.09704,-0.16321,0.09715,-0.24199,0.10886,0.09753,0.02425,-0.022,0.15892,0.20618,-0.23684,-0.14663,0.05369,
			0.11042,1,-0.12822,-0.19143,-0.18822,0.18943,0.09216,-0.31074,-0.18771,-0.0841,0.16604,0.10108,0.10738,-0.02912,0.11395,0.22398,-0.16085,0.13139,-0.15888,0.17561,
			0.15975,-0.12822,1,-0.04698,-0.22569,-0.07349,-0.17824,-0.19955,-0.0821,0.05437,0.02141,0.19926,-0.08384,-0.10899,0.33946,0.0788,0.0612,-0.03674,-0.09496,0.12271,
			0.07215,-0.19143,-0.04698,1,0.05679,-0.05363,-0.03523,0.27162,0.05889,-0.03422,-0.05202,0.24779,0.05109,-0.03347,-0.04865,-0.1377,-0.13698,0.14955,0.08548,-0.02599,
			-0.16194,-0.18822,-0.22569,0.05679,1,-0.19918,-0.01387,0.02249,-0.16517,-0.13269,-0.00939,0.1208,-0.00685,-0.24569,0.13685,-0.03491,-0.23599,-0.11713,0.00331,0.1729,
			0.12573,0.18943,-0.07349,-0.05363,-0.19918,1,0.00275,-0.19135,-0.00738,0.14971,0.02334,0.00178,-0.02642,-0.10979,-0.10028,-0.06943,0.01473,0.00926,0.312,-0.04077,
			0.16026,0.09216,-0.17824,-0.03523,-0.01387,0.00275,1,0.18428,-0.16484,-0.08476,0.03956,-0.19329,-0.17495,0.14818,0.11287,-0.10798,-0.17335,-0.07062,0.1189,-0.13947,
			0.09704,-0.31074,-0.19955,0.27162,0.02249,-0.19135,0.18428,1,-0.20639,0.0125,0.22883,0.09897,-0.18702,0.15289,-0.15701,-0.01498,0.01805,-0.20767,-0.0152,0.01258,
			-0.16321,-0.18771,-0.0821,0.05889,-0.16517,-0.00738,-0.16484,-0.20639,1,-0.21893,0.24765,-0.02113,-0.04481,0.14527,-0.15697,0.18969,0.0092,0.03446,-0.14277,-0.02189,
			0.09715,-0.0841,0.05437,-0.03422,-0.13269,0.14971,-0.08476,0.0125,-0.21893,1,-0.08147,0.05307,0.04591,0.08526,0.16868,0.18528,-0.05258,-0.01021,-0.07407,-0.03105,
			-0.24199,0.16604,0.02141,-0.05202,-0.00939,0.02334,0.03956,0.22883,0.24765,-0.08147,1,0.18587,0.0702,0.096,-0.09226,0.12658,-0.05114,-0.17574,-0.16458,-0.03645,
			0.10886,0.10108,0.19926,0.24779,0.1208,0.00178,-0.19329,0.09897,-0.02113,0.05307,0.18587,1,0.19834,0.09317,0.19145,0.06073,0.05403,0.10893,-0.10967,-0.08809,
			0.09753,0.10738,-0.08384,0.05109,-0.00685,-0.02642,-0.17495,-0.18702,-0.04481,0.04591,0.0702,0.19834,1,-0.12873,0.16739,-0.20469,-0.02687,0.14814,0.19918,-0.11058,
			0.02425,-0.02912,-0.10899,-0.03347,-0.24569,-0.10979,0.14818,0.15289,0.14527,0.08526,0.096,0.09317,-0.12873,1,-0.26324,0.03421,0.02759,-0.13292,0.0148,0.14492,
			-0.022,0.11395,0.33946,-0.04865,0.13685,-0.10028,0.11287,-0.15701,-0.15697,0.16868,-0.09226,0.19145,0.16739,-0.26324,1,0.11387,0.04268,0.03505,0.11692,0.38419,
			0.15892,0.22398,0.0788,-0.1377,-0.03491,-0.06943,-0.10798,-0.01498,0.18969,0.18528,0.12658,0.06073,-0.20469,0.03421,0.11387,1,-0.19631,-0.09071,0.21937,-0.00943,
			0.20618,-0.16085,0.0612,-0.13698,-0.23599,0.01473,-0.17335,0.01805,0.0092,-0.05258,-0.05114,0.05403,-0.02687,0.02759,0.04268,-0.19631,1,0.19802,-0.16292,0.07184,
			-0.23684,0.13139,-0.03674,0.14955,-0.11713,0.00926,-0.07062,-0.20767,0.03446,-0.01021,-0.17574,0.10893,0.14814,-0.13292,0.03505,-0.09071,0.19802,1,-0.05492,-0.1195,
			-0.14663,-0.15888,-0.09496,0.08548,0.00331,0.312,0.1189,-0.0152,-0.14277,-0.07407,-0.16458,-0.10967,0.19918,0.0148,0.11692,0.21937,-0.16292,-0.05492,1,-0.17777,
			0.05369,0.17561,0.12271,-0.02599,0.1729,-0.04077,-0.13947,0.01258,-0.02189,-0.03105,-0.03645,-0.08809,-0.11058,0.14492,0.38419,-0.00943,0.07184,-0.1195,-0.17777,1};
		PetscScalar tmp_b[MATRIX_SIZE] = {   0.26965,
				0.49429,
			   -1.48312,
			   -1.02026,
			   -0.44700,
			    0.10966,
			    1.12874,
			   -0.28996,
			    1.26155,
			    0.47542,
			    1.17412,
			    0.12695,
			   -0.65682,
			   -1.48140,
			    0.15549,
			    0.81855,
			   -0.29259,
			   -0.54079,
			   -0.30864,
			   -1.09659};

		/*the expected solution to system, found with matlab 2014b,
		* error is expected to be below 1e14 for all elements*/
		PetscScalar tmp_x[MATRIX_SIZE] = {73.092608107005603,
										  29.099797385625326,
										 -18.433394494918982,
										 -28.777646223545752,
										   0.478368828160451,
										 -57.825765573842588,
										 -39.926582059511340,
										  -3.547477024603023,
										  35.954640510050552,
										  35.151280818645233,
										  44.057493464710831,
										   8.362334243335667,
										 -67.468535394663505,
										 -17.617047573559240,
										  11.134195015318230,
										 -82.590733482380372,
										 -26.277580602340770,
										  28.729569579045425,
										  81.752960379540212,
										  -2.470184683686883};
										


		this->ke.setSize(MATRIX_SIZE, MATRIX_SIZE);
		this->ke_mapping.allocate(MATRIX_SIZE);

		for(uint32_t ii = 0; ii < MATRIX_SIZE; ++ii) {
			this->rows[ii] = ii;
			this->ke_mapping[ii] = ii;
			this->b[ii] = tmp_b[ii];
			this->x[ii] = tmp_x[ii];
			for(uint32_t jj = 0; jj < MATRIX_SIZE; ++jj){
				this->A[ii*MATRIX_SIZE + jj] = tmp_A[ii*MATRIX_SIZE + jj];
				this->ke(ii, jj) = tmp_A[ii*MATRIX_SIZE + jj];
			}
		}
		/*setup the petsc versions of the marix*/
		PetscErrorCode ierr;
		ierr = VecCreateMPI(PETSC_COMM_WORLD, MATRIX_SIZE, MATRIX_SIZE, &(this->bb));
		EXPECT_EQ(ierr, 0);
		ierr = VecSetOption(this->bb, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
		EXPECT_EQ(ierr, 0);
		ierr = MatCreateSBAIJ(PETSC_COMM_WORLD, 1, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE,
			300, PETSC_NULL, 300, PETSC_NULL, &(this->AA));
		EXPECT_EQ(ierr, 0);
		ierr = MatSetOption(this->AA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		EXPECT_EQ(ierr, 0);
		ierr = MatSetOption(this->AA, MAT_SYMMETRIC, PETSC_TRUE);
		EXPECT_EQ(ierr, 0);
		/*this part will ignore any insertions in the lower triangular, solving
		* all of my problems and making the assembly code 10x less complex */
		ierr = MatSetOption(this->AA, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
		EXPECT_EQ(ierr, 0);
		ierr = KSPCreate(PETSC_COMM_WORLD, &(this->solver));
		EXPECT_EQ(ierr, 0);
		ierr = KSPSetTolerances(this->solver, 1.0e-8, 1e-8,
			PETSC_DEFAULT, 100);
		EXPECT_EQ(ierr, 0);
		ierr = VecDuplicate(this->bb, &(this->xx));
		EXPECT_EQ(ierr, 0);
	}

	virtual void TearDown() {
		MatDestroy(&(this->AA));
		VecDestroy(&(this->xx));
		VecDestroy(&(this->bb));
		KSPDestroy(&(this->solver));
	}

};



TEST_F(MatrixMathTest, SolverTest) {
	PetscErrorCode ierr;

	PetscInt nrows = MATRIX_SIZE;
	MatSetValues(this->AA, nrows, this->rows, nrows, this->rows, this->A, ADD_VALUES);
	VecSetValues(this->bb, nrows, this->rows, this->b, ADD_VALUES);

	ierr = VecAssemblyBegin(this->bb);
	EXPECT_EQ(ierr, 0);
  	ierr = VecAssemblyEnd(this->bb);
  	EXPECT_EQ(ierr, 0);
  	ierr = MatAssemblyBegin(this->AA, MAT_FINAL_ASSEMBLY);
  	EXPECT_EQ(ierr, 0);
  	ierr = MatAssemblyEnd(this->AA, MAT_FINAL_ASSEMBLY);
  	EXPECT_EQ(ierr, 0);

  	/*solve the matrix*/
  	ierr = KSPSetOperators(this->solver, this->AA, this->AA);
  	EXPECT_EQ(ierr, 0);
  	ierr = KSPSetFromOptions(this->solver);
  	EXPECT_EQ(ierr, 0);
  	ierr = KSPSolve(this->solver, this->bb, this->xx);
  	EXPECT_EQ(ierr, 0);


  	PetscScalar out_x[MATRIX_SIZE];

  	VecGetValues(this->xx, nrows, this->rows, out_x);

  	for(uint32_t ii = 0; ii < MATRIX_SIZE; ++ii) {
  		EXPECT_FLOAT_EQ(this->x[ii], out_x[ii]);
  	}
  	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	//MatView(this->AA, PETSC_VIEWER_STDOUT_WORLD);
	//VecView(this->xx, PETSC_VIEWER_STDOUT_WORLD);
}

TEST_F(MatrixMathTest, BulkAssemblerTest) {
	AlgebraicSystem linsys(MATRIX_SIZE);
	std::vector<double> force(MATRIX_SIZE);
	for(uint32_t ii = 0; ii < MATRIX_SIZE; ++ii) {
		force[ii] = this->b[ii];
	}

	linsys.beginAssembly();
	linsys.assemble(this->ke, this->ke_mapping, MATRIX_SIZE);
	linsys.assemble(force, this->ke_mapping, MATRIX_SIZE);
	linsys.synchronize();
	linsys.solve();
	std::vector<double> canidate_x;
	linsys.extractDisplacement(canidate_x);
	for(uint32_t ii = 0; ii < MATRIX_SIZE; ++ii) {
		EXPECT_FLOAT_EQ(this->x[ii] , canidate_x[ii]);
	}
}

TEST_F(MatrixMathTest, AnotherAssemblyInsertionTest) {
	/*here we use the same matrix we have as above, but we 
	* will renumber the degrees of freedom for the global system
	* to simulate the existance of constraints. This way the entire
	* assembly system is tested from end to end, verifying that 
	* all of the mappings between local and global dofs are handled 
	* correctly*/
	uint32_t n_fixed_dofs = 10;
	AlgebraicSystem linsys((MATRIX_SIZE+ n_fixed_dofs));

	/*we choose to add 10 degrees of freedom to the global system that 
	* will be constrained*/
	std::vector<uint64_t> global_constrained_dofs;
	/*the dofs we choose to constrain are arbitrary this test
	* should not matter what dofs are constrained. The only reason
	* we explicitly choose the fixed dofs is so that the test is
	* deterministic*/
	

	global_constrained_dofs.push_back(0); /*make sure minimum is constrained*/
	global_constrained_dofs.push_back(3);
	global_constrained_dofs.push_back(27);	
	global_constrained_dofs.push_back(13);
	global_constrained_dofs.push_back(2);
	global_constrained_dofs.push_back(11);
	global_constrained_dofs.push_back(17);
	global_constrained_dofs.push_back(18);
	global_constrained_dofs.push_back(7);
	global_constrained_dofs.push_back((MATRIX_SIZE + n_fixed_dofs -1)); /*make sure maximum is constrained*/
	ASSERT_EQ(n_fixed_dofs, global_constrained_dofs.size());
	for(uint32_t ii = 0; ii < n_fixed_dofs; ++ii) {
		ASSERT_LT(global_constrained_dofs[ii], (n_fixed_dofs + MATRIX_SIZE));
	}
	/*create associate the fixed dofs with zero displacements*/
	std::vector<double> given_displacements(n_fixed_dofs, 0.0);
	/*relabel the dofs for the matrix, then assemble*/
	linsys.addBoundaryConstraint(given_displacements, global_constrained_dofs);

	linsys.beginAssembly();

	/*relabel the matrix dofs to account for the new fixed members*/
	std::set<PetscInt> taken_dofs;
	for(uint32_t ii = 0; ii < n_fixed_dofs; ++ii) {
		taken_dofs.insert(global_constrained_dofs[ii]);
	}
	/*make sure all of the constrained dofs were unique and 
	* also allow us to test which dofs have already been assigned
	* when relabeling the matrix*/
	ASSERT_EQ(global_constrained_dofs.size(), taken_dofs.size());

	uint32_t curs_indx = 0;
	for(uint32_t ii = 0; ii < (MATRIX_SIZE+n_fixed_dofs); ++ii) {
		if(0 == taken_dofs.count(ii)) {
			/*this dof is not taken*/
			this->ke_mapping[curs_indx] = ii;
			curs_indx++;
		}
	}
	ASSERT_EQ(MATRIX_SIZE, curs_indx);
	/*change the format of the force vector into what method expects*/
	std::vector<double> force(MATRIX_SIZE);
	for(uint32_t ii = 0; ii < MATRIX_SIZE; ++ii) {
		force[ii] = this->b[ii];
	}

	linsys.assemble(this->ke, this->ke_mapping, MATRIX_SIZE);
	linsys.assemble(force, this->ke_mapping, MATRIX_SIZE);
	linsys.synchronize();
	linsys.solve();
	std::vector<double> canidate_x;
	linsys.extractDisplacement(canidate_x);
	uint32_t taken_curs = 0;
	uint32_t not_taken_curs = 0;
	for(uint32_t ii = 0; ii < canidate_x.size(); ++ii) {
		if(1 == taken_dofs.count(ii)) {
			/*make sure that the dofs we fixed have the displacements specified*/
			EXPECT_FLOAT_EQ(given_displacements[taken_curs], canidate_x[ii]);
			taken_curs++;
		} else {
			EXPECT_FLOAT_EQ(this->x[not_taken_curs], canidate_x[ii]);
			not_taken_curs++;
		}
	}
	ASSERT_EQ(MATRIX_SIZE, not_taken_curs);
	ASSERT_EQ(n_fixed_dofs, taken_curs);
}
