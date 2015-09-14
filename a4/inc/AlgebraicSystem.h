#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <map>
#include <stdint.h>

#include <apfDynamicMatrix.h>
#include <canArray.h>
#include <canNewArray.h>

#include <petscksp.h>

#define KNOWN_DOF_MASK ((uint64_t)1<<63)
#define SAMPLE_LOW_32_BITS ((uint32_t)(0xFFFFFFFF))
#define UNMAPPED_VALUE ((uint64_t)-1)

#define SOLVER_ABSOLUTE_TOLERANCE 1.0e-10

class AlgebraicSystem
{
public:
	AlgebraicSystem(std::size_t nGlobalDOFs);
	~AlgebraicSystem();

	void addBoundaryConstraint(
		std::vector<double> const& fixed,
		std::vector<uint64_t> const& node_mapping);
	void beginAssembly();
	void assemble(
		apf::DynamicMatrix const& ke,
		can::NewArray<int> const& node_mapping,
		uint32_t size);
	void assemble(
		std::vector<double> const& fe,
		can::NewArray<int> const& node_mapping,
		uint32_t size);
	PetscErrorCode synchronize();
	PetscErrorCode solve();
	void extractDisplacement(std::vector<double> & disp);
	
	/*these are exposed for unittesting*/
	Mat K;
	Vec F;
	Vec d;
private:
	bool _allow_assembly;
	bool _allow_displacement_extraction;
	/*keep track of which rows have been initialized, this will be
	* used after all assembly is finished to fill zero on the main
	* diagonal of all uninitailized rows which is required by the
	* PETSc solvers we are using */
	std::size_t _num_rows_initialized;
	std::vector<bool> _row_has_been_initialized;
	/*the key is the index in the global numbering of dofs,
	* most significant bit of the value is 1 if the global dof
	* is fixed, with the lower 32 bits being the only valid indices
	* because of interaction with PetscInt not being a reliable 
	* size across platforms*/
	std::map< std::size_t,uint64_t > masks;
	std::vector< double > known_d;
    
    KSP solver;

	std::size_t nGlobalDOFs;
};

#endif
