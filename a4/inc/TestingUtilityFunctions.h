#ifndef TESTING_UTILITY_FUNCTIONS
#define TESTING_UTILITY_FUNCTIONS

#include <gtest/gtest.h>

apf::Vector3 LinearLoad_X(apf::Vector3 const & p)
{
	return apf::Vector3(1000.0, 0, 0);
}

apf::Vector3 LinearLoad_Y(apf::Vector3 const & p) 
{
	return apf::Vector3(0, -1000.0, 0);
}

apf::Vector3 Gravity_Y(apf::Vector3 const & p) 
{
	return apf::Vector3(0, -1000.0, 0);
}


apf::Vector3 SampleLinearBodyForce_Y(apf::Vector3 const & p)
{
	double y_component =  p[1] * 2000.0;
	return apf::Vector3(0, -y_component, 0);
}

enum class MeshTypes { LINEAR_QUAD, QUAD_QUAD, LINEAR_TRI, QUAD_TRI, SERENDIPITY_QUAD };

apf::Mesh2* getMeshFromIndex(MeshTypes index, int X_ELMS, int Y_ELMS, std::string & suffix) {
	MeshBuilder* mesh_builder = new MeshBuilder();
	apf::Mesh2* tmp = NULL;
	switch(index) {
		case MeshTypes::LINEAR_QUAD:
			mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
				0.0, 0.0, 2.0, 2.0);
			suffix = "_linear_quad";
			break;
		case MeshTypes::QUAD_QUAD:
			mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
				0.0, 0.0, 2.0, 2.0);
			apf::changeMeshShape(tmp, apf::getLagrange(2));
			suffix = "_quadratic_quad";
			break;
		case MeshTypes::LINEAR_TRI:
			mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
				0.0, 0.0, 2.0, 2.0);
			suffix = "_linear_tri";
			break;
		case MeshTypes::QUAD_TRI:
			mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
				0.0, 0.0, 2.0, 2.0);
			apf::changeMeshShape(tmp, apf::getLagrange(2));
			suffix = "_quadratic_tri";
			break;
		case MeshTypes::SERENDIPITY_QUAD:
			mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
				0.0, 0.0, 2.0, 2.0);
			apf::changeMeshShape(tmp, apf::getSerendipity());
			suffix = "_quadratic_serendipity";
			break;
		default:
			/*mesh should stay null*/
			break;
	}
	delete mesh_builder;
	return tmp;
}

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

void smallDisplacementNegativeX_2D(
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
		* therefor x terms are the even terms*/
		if(0 == ii %2) {
			dofs[curs_indx] = node_mapping[ii];
			disp[curs_indx] = -1e-2;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}

#endif