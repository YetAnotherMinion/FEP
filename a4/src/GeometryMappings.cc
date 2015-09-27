#include <assert.h>

#include <apfShape.h>

#include "GeometryMappings.h"

void noConstraint(
	apf::MeshEntity *e,
	apf::Mesh *mesh,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> &dofs,
	std::vector<double> & disp)
{
	/*empty the constraint vector*/
	dofs.clear();
	disp.clear();
}

void zeroDisplacementX_2D(
	apf::MeshEntity *e,
	apf::Mesh *mesh,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> & dofs,
	std::vector<double> & disp)
{
	int entity_type = mesh->getType(e);
	uint32_t nnodes = apf::countElementNodes(mesh->getShape(), entity_type);
	/*in 2D we assume there are only 2 dofs per node*/
	apf::NewArray< int > node_mapping(nnodes*2);
	uint32_t tmp_sz = apf::getElementNumbers(nodeNums, e, node_mapping);
	assert((nnodes*2) == tmp_sz);
	/*resize the vector to hold only the fixed dofs in this case is exactly
	* nnodes*/
	dofs.resize(nnodes);
	disp.resize(nnodes);
	std::size_t curs_indx = 0;
	for(std::size_t ii = 0; ii < tmp_sz; ++ii) {
		/*we assume the ordering of local dofs is X_0, Y_0, X_1, Y_1 ...
		* therefor X terms are the even terms*/
		if(0 == ii %2) {
			dofs[curs_indx] = node_mapping[ii];
			disp[curs_indx] = 0.0;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}

void zeroDisplacementY_2D(
	apf::MeshEntity *e,
	apf::Mesh *mesh,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> & dofs,
	std::vector<double> & disp)
{
	int entity_type = mesh->getType(e);
	uint32_t nnodes = apf::countElementNodes(mesh->getShape(), entity_type);
	/*in 2D we assume there are only 2 dofs per node*/
	apf::NewArray< int > node_mapping(nnodes*2);
	uint32_t tmp_sz = apf::getElementNumbers(nodeNums, e, node_mapping);
	assert((nnodes*2) == tmp_sz);
	/*resize the vector to hold only the fixed dofs in this case is exactly
	* nnodes*/
	dofs.resize(nnodes);
	disp.resize(nnodes);
	std::size_t curs_indx = 0;
	for(std::size_t ii = 0; ii < tmp_sz; ++ii) {
		/*we assume the ordering of local dofs is X_0, Y_0, X_1, Y_1 ...
		* therefor y terms are the odd terms*/
		if(1 == ii %2) {
			dofs[curs_indx] = node_mapping[ii];
			disp[curs_indx] = 0.0;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}

GeometryMappings::GeometryMappings()
{
	neumann_map.clear();
	dirchelet_map.clear();
	/*we insert some simple boundary conditions*/
}

GeometryMappings::~GeometryMappings()
{

}

void GeometryMappings::addNeumannMapping(
	uint64_t key,
	apf::Vector3(*fnc_ptr)(apf::Vector3 const& p))
{
	if( 1 == this->neumann_map.count(key)) {
		/*overwrite the existing fnc_ptr if not null,
		* if the function pointer is null, then delete 
		* the key*/
		if(NULL == fnc_ptr) {
			this->neumann_map.erase(key);
		} else {
			this->neumann_map[key] = fnc_ptr;
		}
	} else {
		this->neumann_map.insert(
				std::map<uint64_t, neumann_fnc>::value_type(key, fnc_ptr));
	}
}

void GeometryMappings::addDircheletMapping(
	uint64_t key,
	void(*fnc_ptr)(apf::MeshEntity* e,
					apf::Mesh*,
					apf::Numbering* nodeNums,
					std::vector< uint64_t > & nodes,
					std::vector < double > & d))
{
	if(1 == this->dirchelet_map.count(key)) {
		/*overwrite the existing fnc_ptr if not null,
		* if the function pointer is null, then delete 
		* the key*/
		if(NULL == fnc_ptr) {
			this->dirchelet_map.erase(key);
		} else {
			this->dirchelet_map[key] = fnc_ptr;
		}
	} else {
		this->dirchelet_map.insert(
			std::map<uint64_t, bound_gen>::value_type(key, fnc_ptr));
	}
}

