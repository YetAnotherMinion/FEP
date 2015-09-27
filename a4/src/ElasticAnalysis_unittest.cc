#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>

#include "ElasticAnalysis2D.h"
#include "MeshBuilder.h"
#include "GeometryMappings.h"
#include "AlgebraicSystem.h" /*provides SOLVER_ABSOLUTE_TOLERANCE constant*/

#define YOUNGS_MODULUS  1e8
#define POISSONS_RATIO 0.35

class ElasticAnalysisTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
		mesh_builder = new MeshBuilder();
		mesh = NULL;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		delete mesh_builder;	
	}

};

apf::Vector3 LinearLoad_X(apf::Vector3 const & p)
{
	return apf::Vector3(-1000.0, 0, 0);
}

TEST_F(ElasticAnalysisTest, AppRunTest) {
	mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);
	// mesh_builder->build2DRectTriMesh(this->mesh, 4, 2, 0.0, 0.0, 2.0, 1.0);
	EXPECT_TRUE(this->mesh != NULL);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(this->mesh, apf::getLagrange(2));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	// Nu = POISSONS_RATIO;
	Nu = 0.0;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*Fix the bottom edge in the Y direction only, and the
	* left side in the X direction only*/
	GeometryMappings* geo_map = new GeometryMappings();
	void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
	cnstr_ptr = &zeroDisplacementX_2D;
	geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
	cnstr_ptr = &zeroDisplacementY_2D;
	geo_map->addDircheletMapping(TOP_EDGE, cnstr_ptr);
	//geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

	apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
	traction_ptr = &LinearLoad_X;
	geo_map->addNeumannMapping(RIGHT_EDGE, traction_ptr);

	struct ElasticAnalysisInput input = {
			this->mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	EXPECT_EQ(0, tmp.setup());

	// MatView(tmp.linsys->K, PETSC_VIEWER_STDOUT_WORLD);

	VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);


	EXPECT_EQ(0, tmp.solve());
	// VecView(tmp.linsys->d, PETSC_VIEWER_STDOUT_WORLD);

	EXPECT_EQ(0, tmp.recover());
	std::cout << "=========== Solution ============" << std::endl;
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		std::cout << "d_" << ii << " = " << (tmp.displacement[ii]) << std::endl;
	}
	apf::writeVtkFiles("solution_mesh", this->mesh);
	delete geo_map;
}

TEST_F(ElasticAnalysisTest, PlaneStrainComputation) {
	double E = YOUNGS_MODULUS;
	/*impossible for Nu to go over 0.5*/
	double Nu = POISSONS_RATIO;

	/*manually compute D using different method*/
	double mult_factor = E /((1+Nu)*(1- 2* Nu));

	double elm1 = mult_factor * (1 - Nu);
	double elm2 = mult_factor * (Nu);
	double elm3 = mult_factor * (1 - 2 *Nu) / 2;

	bool use_plane_stress = false;
	apf::Matrix< 3,3 > result_D = buildD(E, Nu, use_plane_stress);

	EXPECT_FLOAT_EQ(elm1, result_D[0][0]);
	EXPECT_FLOAT_EQ(elm2, result_D[0][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[0][2]);
	EXPECT_FLOAT_EQ(elm2, result_D[1][0]);
	EXPECT_FLOAT_EQ(elm1, result_D[1][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[1][2]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][0]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][1]);
	EXPECT_FLOAT_EQ(elm3, result_D[2][2]);
}

TEST_F(ElasticAnalysisTest, PlaneStressComputation) {
	double E = 8e8;
	/*impossible for Nu to go over 0.5*/
	double Nu = 0.35;

	/*manually compute D using different method*/
	double mult_factor = E /(1- Nu * Nu);

	double elm1 = mult_factor ;
	double elm2 = mult_factor * (Nu);
	double elm3 = mult_factor * (1 - Nu) / 2;

	bool use_plane_stress = true;
	apf::Matrix< 3,3 > result_D = buildD(E, Nu, use_plane_stress);

	EXPECT_FLOAT_EQ(elm1, result_D[0][0]);
	EXPECT_FLOAT_EQ(elm2, result_D[0][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[0][2]);
	EXPECT_FLOAT_EQ(elm2, result_D[1][0]);
	EXPECT_FLOAT_EQ(elm1, result_D[1][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[1][2]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][0]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][1]);
	EXPECT_FLOAT_EQ(elm3, result_D[2][2]);

	EXPECT_FLOAT_EQ(result_D[0][1], result_D[1][0]);
	EXPECT_FLOAT_EQ(result_D[2][0], result_D[0][2]);
	EXPECT_FLOAT_EQ(result_D[0][0], result_D[1][1]);
}

class ZeroConstraintZeroTraction : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
		mesh_builder = new MeshBuilder();
		mesh = NULL;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		delete mesh_builder;	
	}

};

TEST_F(ZeroConstraintZeroTraction, LinearQuad) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
	/*now check the displacements*/
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		/*this is a resolution of displacements to the 0.005mm or 5 microns*/
		EXPECT_FLOAT_EQ(0.0, tmp.displacement[ii]);
	}

	delete geo_map;
}

TEST_F(ZeroConstraintZeroTraction, LinearTri) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
	/*now check the displacements*/
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		/*this is a resolution of displacements to the 0.005mm or 5 microns*/
		EXPECT_FLOAT_EQ(0.0, tmp.displacement[ii]);
	}

	delete geo_map;
}

TEST_F(ZeroConstraintZeroTraction, QuadQuad) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	apf::changeMeshShape(mesh, apf::getLagrange(2));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
	/*now check the displacements*/
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		/*this is a resolution of displacements to the 0.005mm or 5 microns*/
		EXPECT_FLOAT_EQ(0.0, tmp.displacement[ii]);
	}

	delete geo_map;
}

TEST_F(ZeroConstraintZeroTraction, SerendipityQuad) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
	/*now check the displacements*/
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		/*this is a resolution of displacements to the 0.005mm or 5 microns*/
		EXPECT_FLOAT_EQ(0.0, tmp.displacement[ii]);
	}

	delete geo_map;
}

TEST_F(ZeroConstraintZeroTraction, QuadTri) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::changeMeshShape(mesh, apf::getLagrange(2));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);
	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
	/*now check the displacements*/
	for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
		EXPECT_FLOAT_EQ(0.0, tmp.displacement[ii]);
	}
	delete geo_map;
}

	// PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	// MatView(linsys.K, PETSC_VIEWER_STDOUT_WORLD);
	// VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);
