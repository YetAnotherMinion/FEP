#include <stdexcept>
#include <climits>
#include <assert.h>
#include <set>

#include <petscmat.h>

#include "AlgebraicSystem.h"

/*when we are bulk copying from petsc to regular std::vectors
* this is the size of the buffer we copy from, this may be a
* very small number, if so it is only for testing,
* IN PRODCUTION change this number so buffer size is close 
* to 1/4 size of cache*/
#define COPY_CHUNK_SIZE 100

AlgebraicSystem::AlgebraicSystem(std::size_t ngd)
	: nGlobalDOFs(ngd)
{
	this->_allow_assembly = false;
	this->_allow_displacement_extraction = false;
	this->_num_rows_initialized = 0;
	/*preallocate the map since its keys will run from [0, nGlobalDOFs-1]*/
	for(std::size_t ii = 0; ii < this->nGlobalDOFs; ++ii) {
		/* -1 is valid intializer for unsigned types since this puts us at
		* UINT_MAX for what ever type we are using. Using this value as 
		* sentinal value for unitialized is valid, since this number
		* represents the indice on a square matrix. Since PETSc,
		* which will be handling this matrix, has never done more than
		* 500 million rows in a matrix before on super computers, it is
		* reasonable to assume that we will never be in danger of overflow
		* for this value for the next 20 years*/
		this->masks.insert(
			std::map< size_t, uint64_t >::value_type(ii, UNMAPPED_VALUE));
	}
	/*make sure there are zero elements in the known_d vector, following
	* code depends on it starting out as zero*/
	this->known_d.clear();
	assert(0 == this->known_d.size());
}

AlgebraicSystem::~AlgebraicSystem()
{
	/*only destroy the matricies if they were allocated*/
	if(true == this->_allow_assembly) {
		MatDestroy(&(this->K));
		VecDestroy(&(this->d));
		VecDestroy(&(this->F));
		KSPDestroy(&(this->solver));
	}
}

void AlgebraicSystem::addBoundaryConstraint(
	std::vector<double> const& fixed,
	std::vector<uint64_t> const& node_mapping)
{
	/*do not add a constraint after assembly has begun*/
	if(this->_allow_assembly == true) {
		throw std::invalid_argument(
			"cannot add constraint after assembly has started");
	}
	std::size_t input_size = fixed.size();
	/*check that there are at least as many mappings as there are elements
	* in fixed that will be mapped. Otherwise we will walk off end of 
	* buffer*/
	if(input_size > node_mapping.size()) {
		throw std::range_error(
			"local force vector rows size exceeded mapping size");
	}
	/*scan over all of the canidate keys for the number of displacements
	* and check if they are already in the mapping. We only check the 
	* displacements given since those are the only ones that actually are
	* added. For example if there are only 3 displacments, but the mapping
	* vector has 4 elements, we have assumed that the first three elements
	* of the mapping correcpond to the first three elements of the mapping
	* In this case if the 4th element of the mapping happens to be a repeated
	* mapping then we actually don't care, since it would have never
	* overwritten anything since there was not a corresponding displacement.
	*/
	std::set<uint64_t> all_possible_keys;
	for(std::size_t ii = 0; ii < input_size; ++ii) {
		uint64_t possible_key = node_mapping[ii];
		/*remember that every valid key was already intialized, so
		* we can also check if this is invalid*/
		if(this->masks.count(possible_key) != 1) {
			throw std::invalid_argument("Global dof mapped is out of range");
		}
		/*TODO FIX: in order to enable the same boundary conditions to be 
		* applied to more than one mesh edge, we need some kind of contention
		* strategy to determine which of the multiple boundary conditions
		* sharing the edge vertex node determine what the dofs endup on that
		* vertex node in question. The easiest thing to do right now is to
		* have every single edge bondary constrain overwrite any exisiting
		* nodes on shared dofs.*/

		// if(this->masks[possible_key] != UNMAPPED_VALUE) {
		// 	throw std::logic_error("Overconstrained error!");
		// }

		/*create a set based on the elements of the node_mapping to check for
		* any repeated keys as well*/
		all_possible_keys.insert(possible_key);
	}
	/*set is used to find duplicate keys*/
	if(all_possible_keys.size() != fixed.size()) {
		throw std::logic_error(
			"Overconstrained inside element mapping, repeated local keys!");
	}
	/*look up current index number of ndogs, we will start
	* inserting each fixed dof at this index growing the vector
	* of known displacements from here*/
	uint64_t cur_ndogs = this->known_d.size();

	for(std::size_t ii = 0; ii < input_size; ++ii) {
		/*TODO FIX: see above, but now it is possible to overwrite mappings
		* so we need to check each key to see if it is already mapped,
		* if so then only change the value of displacement associated instead
		* of pushing a new displacement to the vector */

		/*convert values to fixed width, discarding extra precision
		* if so necessary by using Cs type promotion rules*/
		uint64_t key = node_mapping[ii];
		if(this->masks[key] != UNMAPPED_VALUE) {
			/*find the index in known_d*/
			uint32_t tmp_indx = (this->masks[key] & SAMPLE_LOW_32_BITS);
			std::cout << "key: " << tmp_indx << " already written, !!OVERWRITING" << std::endl;
			this->known_d[tmp_indx] = fixed[ii];

		} else {
			/*look up key an convert it */
			this->masks[key] = cur_ndogs | (KNOWN_DOF_MASK);
			printf("key %lx \n", this->masks[key]);
			++cur_ndogs;
			/*add the proscribed displacements to our known_d(isplacement) store*/
			this->known_d.push_back(fixed[ii]);
		}
	}

}

void AlgebraicSystem::beginAssembly() {
	/*allow beginAssembly to be called only once per object
	* lifetime*/
	if(this->_allow_assembly) {
		throw std::logic_error(
			"Inconsistant state, cannot call beginAssembly more than once");
	}
	_allow_assembly = true;
	/*now compute the remaining degrees of freedom*/
	std::size_t ndogs = this->known_d.size();
	std::size_t n_eq = this->nGlobalDOFs - ndogs;
	/*initialize the flag of */
	this->_row_has_been_initialized.clear();
	this->_row_has_been_initialized.resize(n_eq, false);
	/*go throught the map and map the free degrees of freedom
	* to indicies that will be used with Petsc matricies and vectors
	* to solve the remaining degrees of freedom*/
	uint64_t n_eq_counter = 0;
	for(std::size_t ii = 0; ii < this->nGlobalDOFs; ++ii) {
		if(this->masks[ii] == UNMAPPED_VALUE) {
			/*clear the most significant bit to indicate that this is
			* a free dof */
			this->masks[ii] = n_eq_counter & (~(KNOWN_DOF_MASK));
			++n_eq_counter;
		} else {
			/*this dof was already labled*/
		}
	}
	/*we should have seen exactly n_eq free dofs*/
	assert(n_eq == n_eq_counter);

	/*we construct the stiffness matrix only for degress of freedom
	* which have not been fixed, 
	* WARNING: this program is intended to run on a single process
	* there for the local size of all petsc variables is set to
	* the same as the global size.*/
	PetscErrorCode ierr;
	ierr = VecCreateMPI(PETSC_COMM_WORLD, n_eq, n_eq, &(this->F));
	assert(0 == ierr);
	ierr = VecSetOption(this->F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
	assert(0 == ierr);
	ierr = MatCreateSBAIJ(PETSC_COMM_WORLD, 1, n_eq, n_eq, n_eq, n_eq,
		300, PETSC_NULL, 300, PETSC_NULL, &(this->K));
	assert(0 == ierr);
	ierr = MatSetOption(this->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	assert(0 == ierr);
	ierr = MatSetOption(this->K, MAT_SYMMETRIC, PETSC_TRUE);
	assert(0 == ierr);
	/*this part will ignore any insertions in the lower triangular, solving
	* all of my problems and making the assembly code 10x less complex */
	ierr = MatSetOption(this->K, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
	assert(0 == ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD, &(this->solver));
	assert(0 == ierr);
	ierr = KSPSetTolerances(this->solver, 1.0e-10, SOLVER_ABSOLUTE_TOLERANCE,
		PETSC_DEFAULT, 100);
	assert(0 == ierr);
	ierr = VecDuplicate(this->F, &(this->d));
	assert(0 == ierr);
}

void AlgebraicSystem::assemble(
	apf::DynamicMatrix const& ke,
	can::NewArray<int> const& node_mapping,
	uint32_t nLocalDOFs)
{
	/*CHECK if assembly is possible*/
	if(false == this->_allow_assembly) {
		throw std::invalid_argument("must call beginAssembly before assembly");
	}

	/*iterate over the input matrix and add each value
	* into the global matrix according to the mapping in */
	std::size_t nrows = ke.getRows();
	std::size_t ncol = ke.getColumns();
	/*sanity checks that we actually have a mapglobal_n_dofsping*/
	if(nLocalDOFs > this->nGlobalDOFs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*this will accept a non square matrix as long as its mappings fit*/
	/*check that every row will have a mapping*/
	if((nrows > this->nGlobalDOFs) || (nrows > nLocalDOFs)) {
		throw std::range_error("local rows size exceeded mapping size");
	}
	/*check that every column will have mapping and will fit into the
	* global matrix*/
	if((ncol > this->nGlobalDOFs) || (ncol > nLocalDOFs)) {
		throw std::range_error("local column size exceeded mapping size");
	}
	/*it will be cheaper to allocate the full possible size for the
	* global indice vector than to compute exact size */
	PetscInt *idxM = new PetscInt[nrows];
	PetscInt *idxN = new PetscInt[ncol];
	PetscInt *ix = new PetscInt[nrows];

	/*determine size of the assembly matrix & any force vector contributions*/
	
	/*cache the local indicies we will need*/
	std::vector< std::size_t > keep_rows, keep_cols, add_to_f_col;
	keep_rows.clear();
	keep_cols.clear();
	add_to_f_col.clear();
	std::vector< double > local_displacements;
	local_displacements.clear();
	std::vector<uint32_t> rows_marked_valid;
	rows_marked_valid.clear();

	/*find the indices of the element stiffness we will keep*/
	PetscInt idxM_size = 0;

	for(std::size_t ii = 0; ii < nrows; ++ii){
		/*get the reduced global index using the map*/
		int tmp_global_index = node_mapping[ii];
		if(tmp_global_index < 0){
			throw std::logic_error(
				"cannot have negative indices for global dofs");
		}
		/*convert to unsigned type using promotion*/
		std::size_t tmp_key = tmp_global_index;
		uint64_t tmp_val = this->masks[tmp_key];
		/*test if dof is fixed or not*/
		if(!(tmp_val & KNOWN_DOF_MASK)){
			/*extract the 32 bit number*/
			uint32_t tmp_mapped_val = (tmp_val & SAMPLE_LOW_32_BITS);
			idxM[idxM_size] = tmp_mapped_val;
			/*mark this row as initialized, we must be careful with this
			* because we must also make sure that at least one column is
			* a contributor to the stiffness matrix, otherwise we could
			* only be contributing to the force vector. */
			if(this->_row_has_been_initialized[tmp_mapped_val]) {
				/*do nothing if it has already been passed*/
			} else {
				this->_row_has_been_initialized[tmp_mapped_val] = true;
				this->_num_rows_initialized++;
				/*we save each row we marked as valid, so we can roll back
				* any marked rows if it turns out there are no contributions
				* to the stiffness matrix for this element*/
				rows_marked_valid.push_back(tmp_mapped_val);
			}
			
			//std::cout << "ii: " << ii << " maps to: " << idxM[idxM_size] << std::endl;
			++idxM_size;
			keep_rows.push_back(ii);
		} else {
			//std::cout << "dof already mapped" << std::endl;
		}
	}

	/*find the columns of the element stiffness matrix we will keep and 
	* the columns that will contribute to the force vector*/
	PetscInt idxN_size = 0;
	PetscInt ni = 0;
	for(std::size_t ii = 0; ii < ncol; ++ii) {
		/*get the reduced global index using the map*/
		int tmp_global_index = node_mapping[ii];
		if(tmp_global_index < 0){
			throw std::logic_error(
				"cannot have negative indices for global dofs");
		}
		/*convert to unsigned type using promotion*/
		std::size_t tmp_key = tmp_global_index;
		uint64_t tmp_val = this->masks[tmp_key];
		/*test if dof is fixed or not*/
		if(!(tmp_val & KNOWN_DOF_MASK)){
			/*high bit not set means it is a free degree*/
			/*extract the 32 bit number*/
			idxN[idxN_size] = (PetscInt) (tmp_val & SAMPLE_LOW_32_BITS);
			//std::cout << "ii: " << ii << " maps to: " << idxM[idxM_size] << std::endl;
			++idxN_size;
			keep_cols.push_back(ii);
		} else {
			/*high bit IS set, so it is fixed*/

			ix[ni] = (PetscInt) (tmp_val & SAMPLE_LOW_32_BITS);
			add_to_f_col.push_back(ii);
			local_displacements.push_back(this->known_d[(tmp_val & SAMPLE_LOW_32_BITS)]);
			
			std::cout << "col contributes to force vector: " << ii << " mapping: " << node_mapping[ii] <<  std::endl;
			
			++ni;
		}
	}
	/*each of the PetscInt's will be larger than zero*/
	assert(static_cast<uint64_t>(idxM_size) == keep_rows.size());
	assert(static_cast<uint64_t>(idxN_size) == keep_cols.size());
	assert(static_cast<uint64_t>(ni) == local_displacements.size());
	assert(static_cast<uint64_t>(ni) == add_to_f_col.size());
	/*here is where we roll back any marked valid rows if in fact
	* non contribute, aka when keep_cols.size == 0*/
	if(0 == keep_cols.size()){
		for(std::size_t ii = 0; ii < rows_marked_valid.size(); ++ii) {
			this->_row_has_been_initialized[rows_marked_valid[ii]] = false;
			this->_num_rows_initialized--;
		}
	}

	/*allocate a two dimensional storage for stiffness matrix 
	* using PetscInts which by the default options chosen elsewhere
	* should be row-oriented*/
	PetscScalar *values = new PetscScalar[idxM_size * idxN_size];
	PetscScalar *y = new PetscScalar[idxM_size];

	/*loop over element stiffness matricies and copy contents into 
	* the petsc approved formats*/

	for(std::size_t ii = 0; ii < keep_rows.size(); ++ii){
		/*skip all rows that are already determined*/
		for(std::size_t jj = 0; jj < keep_cols.size(); ++jj){
			/*get the indicies we want in the local stiffness matrix*/
			uint64_t kk = keep_rows[ii];
			uint64_t ll = keep_cols[jj];
			/*load the petsc scalar values which are known to be not zero
			* by virtue of check above from the double*/
			PetscScalar tmp = ke(kk,ll);
			std::size_t computed_index = ii * idxM_size + jj;
			values[computed_index] = tmp;
		}
		/*move the stiffness contribution for the fixed dofs to the force
		* vector side*/
		for(std::size_t jj = 0; jj < add_to_f_col.size(); ++jj) {
			uint64_t kk = keep_rows[ii];
			uint64_t ll = add_to_f_col[jj];
			/*load the petsc scalar values which are known to be not zero
			* by virtue of check above from the double*/
			PetscScalar tmp = ke(kk,ll);
			tmp *= ( -1 * local_displacements[jj]);
			y[jj] = tmp;
		}
	}
	/*now that all of the above has been computed*/
	MatSetValues(this->K, idxM_size, idxM, idxN_size, idxN, values, ADD_VALUES);
	VecSetValues(this->F, ni, ix, y, ADD_VALUES);

	delete[] y;
	delete[] ix;
	delete[] idxM;
	delete[] idxN;
	delete[] values;
}

void AlgebraicSystem::assemble(
	std::vector<double> const& fe,
	can::NewArray<int> const& node_mapping,
	uint32_t nLocalDOFs)
{
	/*iterate over the input vector and add each value
	* into the global force vector according to the mapping in */
	std::size_t nrows = fe.size();
	/*sanity checks that we actually have a mapglobal_n_dofsping*/
	if(nLocalDOFs > this->nGlobalDOFs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*check that every row will have a mapping*/
	if((nrows > this->nGlobalDOFs) || (nrows > nLocalDOFs)) {
		throw std::range_error("local rows size exceeded mapping size");
	}
	std::size_t f_size = fe.size();
	std::cout << "f_size" << f_size << std::endl;

	PetscInt *ix = new PetscInt[f_size];
	PetscScalar *y = new PetscScalar[f_size];
	PetscInt ni = 0;

	for(std::size_t ii = 0; ii < f_size; ++ii) {
		/*get the reduced global index using the map*/
		int tmp_global_index = node_mapping[ii];
		if(tmp_global_index < 0){
			throw std::logic_error(
				"cannot have negative indices for global dofs");
		}
		/*convert to unsigned type using promotion*/
		std::size_t tmp_key = tmp_global_index;
		uint64_t tmp_val = this->masks[tmp_key];
		/*test if dof is fixed or not*/
		if(!(tmp_val & KNOWN_DOF_MASK)){
			/*high bit not set means it is a free degree*/
			/*extract the 32 bit number*/
			ix[ni] = static_cast<PetscInt>((tmp_val & SAMPLE_LOW_32_BITS));
			PetscScalar tmp = fe[ii];
			y[ii] = tmp;
			++ni;
		} else {
			std::cout << "dropping force contribution" << std::endl;
			/*high bit IS set, so it is fixed, so the force contribution is
			* silently dropped */
		}
	}
	VecSetValues(this->F, ni, ix, y, ADD_VALUES);

	delete[] ix;
	delete[] y;
}

PetscErrorCode AlgebraicSystem::synchronize()
{
	/*CHECK if assembly is possible*/
	if(false == this->_allow_assembly) {
		throw std::invalid_argument("must call beginAssembly before synchronize");
	}
	PetscErrorCode ierr;
	/*allocate a buffer of rows that the main diagonal
	* must be initialized to zero*/
	PetscInt numRows = this->_row_has_been_initialized.size()
						- this->_num_rows_initialized;
	PetscInt rows_to_zero; 
	PetscScalar values;
	std::size_t curs_indx = 0;
	for(std::size_t ii = 0; ii < this->_row_has_been_initialized.size(); ++ii) {
		if(false == this->_row_has_been_initialized[ii]) {
			rows_to_zero = ii;
			values = 0.0;
			/*TODO: FIX: right now just set the value to zero of the diagonal 
			* element, so that other tests pass because they do not set the 
			* the other rows, and removing all of the rows which are not set messes
			* these tests up bu shrinking matrix size to zero*/
			MatSetValues(this->K, 1, &rows_to_zero, 1, &rows_to_zero,
				&values, ADD_VALUES);
			// ierr = MatZeroRows(this->K, numRows, rows_to_zero, 0.0, this->d, this->F);
			curs_indx++;
		}
	}
	assert(curs_indx == static_cast<uint64_t>(numRows));

	ierr = VecAssemblyBegin(this->F);
  	CHKERRQ(ierr);
  	ierr = VecAssemblyEnd(this->F);
  	CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
  	CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
  	return (PetscErrorCode)0;
}

PetscErrorCode AlgebraicSystem::solve()
{
	/*CHECK if assembly is possible*/
	if(false == this->_allow_assembly) {
		throw std::invalid_argument("must call beginAssembly before solving");
	}
	std::cout << "Solving system of: "
		<< this->masks.size() - this->known_d.size() 
		<< " dofs" << std::endl;
	PetscErrorCode ierr;
	/*we use same matrix as preconditioner*/
	ierr = KSPSetOperators(this->solver, this->K, this->K);
	CHKERRQ(ierr);
  	ierr = KSPSetFromOptions(this->solver);
  	CHKERRQ(ierr);
  	ierr = KSPSolve(this->solver, this->F, this->d);
  	CHKERRQ(ierr);
  	/*zero is no error*/
  	this->_allow_displacement_extraction = true;
  	return (PetscErrorCode)0;
}

void AlgebraicSystem::extractDisplacement(std::vector<double> & disp)
{
	if(false == this->_allow_displacement_extraction) {
		throw std::invalid_argument(
			"must call Solve before extracting displacements");
	}
	disp.resize(this->nGlobalDOFs, 0.0);
	/*get values from vector in chunks*/

	PetscScalar *y = new PetscScalar[COPY_CHUNK_SIZE];
	PetscInt *ix = new PetscInt[COPY_CHUNK_SIZE];

	/*iterate over the map, so we can fill in relevant known displacements
	* as we get to them without having to break up the vector reading 
	* stream into ireguarly sized chunks */
	std::map< std::size_t,uint64_t >::iterator map_it = this->masks.begin();

	for(std::size_t ii = 0; ii < this->nGlobalDOFs;) {
		if(this->masks.end() == map_it) {
			throw std::logic_error("unmapping went horribly wrong");
		}
		/*save a copy of where the map iterator was to fill it in*/
		std::map< std::size_t,uint64_t >::iterator back_fill = map_it;
		/*we do not assume that the map iterator is in order of keys*/
		/*keep this counter globally visible so we can reuse it when filling
		* in the output vector*/
		std::size_t jj;
		for(jj = 0; jj < COPY_CHUNK_SIZE;) {
			/*attempt to find the next free degree, if the degree
			* is actually fixed then just add it to correct place
			* in output vector and move on to the next key in the
			* map*/
			uint64_t tmp_val = map_it->second;
			/*test if dof is fixed or not*/
			if(!(tmp_val & KNOWN_DOF_MASK)){
				/*high bit not set means it is a free degree*/
				/*extract the 32 bit number*/
				ix[jj] = static_cast<PetscInt>((tmp_val & SAMPLE_LOW_32_BITS));
				/*because a dof was found, add its index to be retrieved
				* from the petsc vector*/
				++jj;
			} else {
			}
			/*increment to next key*/
			map_it++;
			if(map_it == this->masks.end()) {
				/*when we run out of mappings we stop immediatly*/
				break;
			}
		}
		/*jj is now the number of free dofs transfered*/
		VecGetValues(this->d, static_cast<PetscInt>(jj), ix, y);
		/*now using the back fill iterator fill in the remaining values*/
		for(std::size_t kk = 0; back_fill != map_it;) {
			std::size_t tmp_key = back_fill->first;
			uint64_t tmp_val = back_fill->second;
			/*test if dof is fixed or not*/
			if(!(tmp_val & KNOWN_DOF_MASK)){
				/*high bit not set means it is a free degree*/
				assert(kk < jj);
				disp[tmp_key] = y[kk];
				/*move on to the next element in the buffer retrived
				* from the buffer vector*/
				++kk;
			} else {
				disp[tmp_key] = this->known_d[(tmp_val & SAMPLE_LOW_32_BITS)];
			}
			/*indicate that we wrote one dof*/
			++ii;
			/*increment to next key*/
			back_fill++;
		}
	}
	delete[] ix;
	delete[] y;
}
