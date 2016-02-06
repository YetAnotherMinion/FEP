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
	double y_component =  p[1] * 1000;
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

#endif