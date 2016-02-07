#include <assert.h>

#include <apfShape.h>

#include "GeometryMappings.h"


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

