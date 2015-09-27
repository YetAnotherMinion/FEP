#ifndef GEOMETRY_MAPPING_H
#define GEOMETRY_MAPPING_H 1

#include <map>
#include <vector>
#include <set>
#include <stdint.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>

#define VERT_BC_TAG_NAME "vertBCtag"
#define EDGE_BC_TAG_NAME "edgeBCtag"
#define FACE_BC_TAG_NAME "faceBCtag"


/*include some simple defaults*/
enum GeoMapOptions {
	NONE = 0,
	ZERO_DISP_X = 1,
	ZERO_DISP_Y = 2,
	FIXED_XY = 3,
	GRAVITY_BODY_FORCE = 4,
	TRACTION_1 = 5,
	TRACTION_2 = 6,
	TRACTION_3 = 7,
	ALL_FACES = 99,
	LEFT_EDGE = 100, /*these are for 2D rectangle meshes*/
	RIGHT_EDGE = 101,
	TOP_EDGE = 102,
	BOT_EDGE = 103,
	VERT_X0Y0 = 104,
	VERT_X1Y0 = 105,
	VERT_X0Y1 = 106,
	VERT_X1Y1 = 107
};

/*add a default zero constraint option*/
void noConstraint(
	apf::MeshEntity*,
	apf::Mesh*,
	apf::Numbering*,
	std::vector<uint64_t> &,
	std::vector<double> &);
void zeroDisplacementX_2D(
	apf::MeshEntity*,
	apf::Mesh*,
	apf::Numbering*,
	std::vector<uint64_t> &,
	std::vector<double> &);
void zeroDisplacementY_2D(
	apf::MeshEntity*,
	apf::Mesh*,
	apf::Numbering*,
	std::vector<uint64_t> &,
	std::vector<double> &);


typedef apf::Vector3(*neumann_fnc)(apf::Vector3 const& );
typedef void(*bound_gen)(apf::MeshEntity*, apf::Mesh*, apf::Numbering* nodeNums, std::vector< uint64_t >&, std::vector < double > & );

class GeometryMappings
{
	/*manually configured in the .cc file, this applies
	* all boundary conditions to geometric model based on integer tags
	*/

public:
	GeometryMappings(apf::Mesh *m);
	~GeometryMappings();

	void addNeumannMapping(uint64_t key, neumann_fnc);
	void addDircheletMapping(uint64_t key, bound_gen);

	std::map< uint64_t,neumann_fnc > neumann_map;
	std::map< uint64_t,bound_gen > dirchelet_map;

private:
	apf::Mesh* mesh;

};


#endif

