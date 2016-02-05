#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <string>

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


TEST_F(ElasticAnalysisTest, AppRunTest) {
	mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);
	// mesh_builder->build2DRectTriMesh(this->mesh, 4, 2, 0.0, 0.0, 2.0, 1.0);
	EXPECT_TRUE(this->mesh != NULL);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(this->mesh, apf::getLagrange(2));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	//Nu = 0.0;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*Fix the bottom edge in the Y direction only, and the
	* left side in the X direction only*/
	GeometryMappings* geo_map = new GeometryMappings();
	void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
	cnstr_ptr = &zeroDisplacementX_2D;
	geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
	cnstr_ptr = &zeroDisplacementY_2D;
	//geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
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

	// VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);


	EXPECT_EQ(0, tmp.solve());
	// VecView(tmp.linsys->d, PETSC_VIEWER_STDOUT_WORLD);

	EXPECT_EQ(0, tmp.recover());
	// std::cout << "=========== Solution ============" << std::endl;
	// for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
	// 	std::cout << "d_" << ii << " = " << (tmp.displacement[ii]) << std::endl;
	// }
	
	// apf::writeASCIIVtkFiles("solution_mesh", this->mesh);
	delete geo_map;
}

TEST_F(ElasticAnalysisTest, SideConcaveTest) {
	mesh_builder->build2DRectQuadMesh(this->mesh, 6, 6, 0.0, 0.0, 2.0, 1.0);
	// mesh_builder->build2DRectTriMesh(this->mesh, 4, 2, 0.0, 0.0, 2.0, 1.0);
	EXPECT_TRUE(this->mesh != NULL);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(this->mesh, apf::getLagrange(2));
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	//Nu = 0.0;
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
	geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

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
	// VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);


	EXPECT_EQ(0, tmp.solve());
	// VecView(tmp.linsys->d, PETSC_VIEWER_STDOUT_WORLD);

	EXPECT_EQ(0, tmp.recover());

	// std::cout << "=========== Solution ============" << std::endl;
	// for(std::size_t ii = 0; ii < tmp.displacement.size(); ++ii) {
	// 	std::cout << "d_" << ii << " = " << (tmp.displacement[ii]) << std::endl;
	// }
	
	apf::writeASCIIVtkFiles("side_concave", this->mesh);
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
	double Nu = 0.35;
	/*impossible for Nu to go over 0.5*/
	ASSERT_LE(Nu, 0.5);

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

enum class MeshTypes { LINEAR_QUAD, QUAD_QUAD, LINEAR_TRI, QUAD_TRI, SERENDIPITY_QUAD };

class SampleProblems : public ::testing::Test,
	public ::testing::WithParamInterface<MeshTypes> 
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;
	std::string suffix;

	apf::Mesh2* getMeshFromIndex(MeshTypes index, int X_ELMS, int Y_ELMS) {
		apf::Mesh2* tmp = NULL;
		switch(index) {
			case MeshTypes::LINEAR_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				this->suffix = "_linear_quad";
				break;
			case MeshTypes::QUAD_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getLagrange(2));
				this->suffix = "_quadratic_quad";
				break;
			case MeshTypes::LINEAR_TRI:
				this->mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				this->suffix = "_linear_tri";
				break;
			case MeshTypes::QUAD_TRI:
				this->mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getLagrange(2));
				this->suffix = "_quadratic_tri";
				break;
			case MeshTypes::SERENDIPITY_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getSerendipity());
				this->suffix = "_quadratic_serendipity";
				break;
			default:
				/*mesh should stay null*/
				break;
		}
		return tmp;
	} 

	virtual void SetUp() {
		/*use mesh builder to abstract mesh creations*/
		mesh_builder = new MeshBuilder();
		mesh = NULL;
		suffix = "";
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
		delete mesh_builder;
	}
};

INSTANTIATE_TEST_CASE_P(GeometryTags, SampleProblems,
	::testing::Values(MeshTypes::LINEAR_QUAD, MeshTypes::QUAD_QUAD,
					  MeshTypes::LINEAR_TRI, MeshTypes::QUAD_TRI,
					  MeshTypes::SERENDIPITY_QUAD));

TEST_P(SampleProblems, Gravity) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 80;
	uint32_t Y_ELMS = 80;
	this->mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	//Nu = 0.0;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*Fix the bottom edge in the Y direction only, and the
	* left side in the X direction only*/
	GeometryMappings* geo_map = new GeometryMappings();
	void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
	cnstr_ptr = &zeroDisplacementX_2D;
	geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
	cnstr_ptr = &zeroDisplacementY_2D;
	geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

	apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
	traction_ptr = &Gravity_Y;
	geo_map->addNeumannMapping(ALL_FACES, traction_ptr);

	struct ElasticAnalysisInput input = {
			this->mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	// PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	EXPECT_EQ(0, tmp.setup());

	// MatView(tmp.linsys->K, PETSC_VIEWER_STDOUT_WORLD);

	// VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);


	EXPECT_EQ(0, tmp.solve());
	// VecView(tmp.linsys->d, PETSC_VIEWER_STDOUT_WORLD);
	EXPECT_EQ(0, tmp.recover());
	std::cout << "************************************" << std::endl;
	std::cout << "*\t" << this->suffix << std::endl;
	std::cout << "************************************" << std::endl;
	std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	std::cout << std::endl;
	std::string out_name = "gravity";
	out_name += this->suffix;
	
	// apf::writeASCIIVtkFiles(out_name.c_str(), this->mesh);
	apf::writeVtkFiles(out_name.c_str(), this->mesh);
	delete geo_map;
}

TEST_P(SampleProblems, ZeroConstraintZeroTraction) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 4;
	uint32_t Y_ELMS = 3;
	this->mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
	ASSERT_TRUE(this->mesh != NULL);
	/*physical parameters*/
	double E, Nu;
	E = YOUNGS_MODULUS;
	Nu = POISSONS_RATIO;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			this->mesh,
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

	for(auto pair : tmp.strain) {
		/*each strain component should be zero*/
		for(uint32_t ii; ii < 3; ++ii) {
			EXPECT_FLOAT_EQ(0.0, pair.second[ii]);
		}
		// std::cout << "x = " << pair.first << ":" << std::endl;
		// std::cout << "\te = " << pair.second << std::endl << std::endl;
	};
	EXPECT_FLOAT_EQ(0.0, tmp.strain_energy);

	std::cout << "************************************" << std::endl;
	std::cout << "*\t" << this->suffix << std::endl;
	std::cout << "************************************" << std::endl;
	std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	std::cout << std::endl;
	std::string out_name = "zero_zero";
	out_name += this->suffix;
	
	apf::writeASCIIVtkFiles(out_name.c_str(), this->mesh);
	// apf::writeVtkFiles(out_name.c_str(), this->mesh);
	delete geo_map;
}

#define LINEAR_X_LOAD 1000.0

apf::Vector3 SampleLinearLoad_X(apf::Vector3 const & p)
{
	return apf::Vector3(LINEAR_X_LOAD, 0, 0);
}

apf::Vector3 SampleLinearBodyForce_Y(apf::Vector3 const & p)
{
	double y_component = 500.0 + p[0] * 50;
	return apf::Vector3(0, -y_component, 0);
}

TEST_P(SampleProblems, LinearTraction) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 2;
	uint32_t Y_ELMS = 1;
	this->mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
	ASSERT_TRUE(this->mesh != NULL);
	/*physical parameters*/
	double E, Nu;
	E = 8e8;
	Nu = 0.35;
	uint32_t integration_order = 4;
	bool reorder_flag = true;

	GeometryMappings* geo_map = new GeometryMappings();

	void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
	cnstr_ptr = &zeroDisplacementX_2D;
	geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
	cnstr_ptr = &zeroDisplacementY_2D;
	//geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
	//geo_map->addDircheletMapping(TOP_EDGE, cnstr_ptr);
	geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

	apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
	traction_ptr = &SampleLinearLoad_X;
	geo_map->addNeumannMapping(RIGHT_EDGE, traction_ptr);

	struct ElasticAnalysisInput input = {
			this->mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());

	/*expect uniform stresses*/
	for(auto pair : tmp.strain) {
		/*linear test should have strain in perpendicular direction
		* following Hook's law*/
		EXPECT_NEAR(pair.second[0] * -Nu, pair.second[1], 1e-13);
		/* zero strain in Z to satisfy our intial assumptions*/
		EXPECT_NEAR(0.0, pair.second[2], 1e-13);
	}
	for(auto pair : tmp.stress) {
		EXPECT_FLOAT_EQ(LINEAR_X_LOAD, pair.second[0]);
	}

	// std::cout << "************************************" << std::endl;
	// std::cout << "*\t\t" << this->suffix << std::endl;
	// std::cout << "************************************" << std::endl;
	// std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	// std::cout << std::endl;
	// std::string bar = "linear_load";
	// bar += this->suffix;

	delete geo_map;
}
#define UPPER_H_LIMIT 12

TEST_P(SampleProblems, StrainEnergyConstantBodyForce) {
	MeshTypes index = GetParam();

	double energies[UPPER_H_LIMIT];
	int problem_dofs[UPPER_H_LIMIT];

	for(uint32_t h = 1; h < UPPER_H_LIMIT; ++h) {
		uint32_t X_ELMS = h;
		uint32_t Y_ELMS = h+1;
		apf::Mesh2* mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
		std::cout << "\t\t\t\t" << this->suffix << std::endl;
		ASSERT_TRUE(mesh != NULL);
		/*physical parameters*/
		double E, Nu;
		E = 8e8;
		Nu = 0.35;
		uint32_t integration_order = 4;
		bool reorder_flag = true;

		GeometryMappings* geo_map = new GeometryMappings();

		void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
		cnstr_ptr = &zeroDisplacementX_2D;
		geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
		cnstr_ptr = &zeroDisplacementY_2D;
		//geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
		//geo_map->addDircheletMapping(TOP_EDGE, cnstr_ptr);
		geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

		apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
		traction_ptr = &Gravity_Y;
		geo_map->addNeumannMapping(ALL_FACES, traction_ptr);
		// traction_ptr = &SampleLinearLoad_X;
		// geo_map->addNeumannMapping(RIGHT_EDGE, traction_ptr);

		struct ElasticAnalysisInput input = {
				mesh,
				geo_map,
				integration_order,
				E,
				Nu,
				reorder_flag};

		ElasticAnalysis2D *tmp = new ElasticAnalysis2D(input);

		EXPECT_EQ(0, tmp->setup());
		EXPECT_EQ(0, tmp->solve());
		EXPECT_EQ(0, tmp->recover());

		problem_dofs[h] = tmp->linsys->freeDOFs;
		energies[h] = tmp->strain_energy;
		/*preserve a record of some of these meshes*/
		if( 10 == h ) {
			std::string out_name("constant_body_force_2x2");
			out_name += this->suffix;
			apf::writeASCIIVtkFiles(out_name.c_str(), mesh);
		}

		delete tmp;
		mesh->destroyNative();
		apf::destroyMesh(mesh);
		delete geo_map;
	}
	std::cout << "constant_body_force_2x2" << this->suffix << " = [" << std::endl; 
	for(auto ii = 1; ii < UPPER_H_LIMIT; ++ii) {
		std::cout << ii << ", " << problem_dofs[ii] << ", ";
		std::cout << std::setprecision(20);
		std::cout << energies[ii] <<  std::setprecision(5) << std::endl;
	}
	std::cout << "];" << std::endl;
}

TEST_P(SampleProblems, StrainEnergyLinearBodyForce) {
	MeshTypes index = GetParam();

	double energies[UPPER_H_LIMIT];
	int problem_dofs[UPPER_H_LIMIT];

	for(uint32_t h = 1; h < UPPER_H_LIMIT; ++h) {
		uint32_t X_ELMS = h;
		uint32_t Y_ELMS = h+1;
		apf::Mesh2* mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
		std::cout << "\t\t\t\t" << this->suffix << std::endl;
		ASSERT_TRUE(mesh != NULL);
		/*physical parameters*/
		double E, Nu;
		E = 8e8;
		Nu = 0.35;
		uint32_t integration_order = 4;
		bool reorder_flag = true;

		GeometryMappings* geo_map = new GeometryMappings();

		void (*cnstr_ptr)(apf::MeshEntity*, apf::Mesh*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
		cnstr_ptr = &zeroDisplacementX_2D;
		geo_map->addDircheletMapping(LEFT_EDGE, cnstr_ptr);
		cnstr_ptr = &zeroDisplacementY_2D;
		geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

		apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
		traction_ptr = &SampleLinearBodyForce_Y;
		geo_map->addNeumannMapping(ALL_FACES, traction_ptr);

		struct ElasticAnalysisInput input = {
				mesh,
				geo_map,
				integration_order,
				E,
				Nu,
				reorder_flag};

		ElasticAnalysis2D *tmp = new ElasticAnalysis2D(input);

		EXPECT_EQ(0, tmp->setup());
		EXPECT_EQ(0, tmp->solve());
		EXPECT_EQ(0, tmp->recover());

		problem_dofs[h] = tmp->linsys->freeDOFs;
		energies[h] = tmp->strain_energy;
		/*preserve a record of some of these meshes*/
		if( 10 == h ) {
			std::string out_name("linear_body_force_2x2");
			out_name += this->suffix;
			apf::writeASCIIVtkFiles(out_name.c_str(), mesh);
		}


		delete tmp;
		mesh->destroyNative();
		apf::destroyMesh(mesh);
		delete geo_map;
	}

	std::cout << "linear_body_force" << this->suffix << " = [" << std::endl; 
	for(auto ii = 1; ii < UPPER_H_LIMIT; ++ii) {
		std::cout << ii << ", " << problem_dofs[ii] << ", ";
		std::cout << std::setprecision(20);
		std::cout << energies[ii] <<  std::setprecision(5) << std::endl;
	}
	std::cout << "];" << std::endl;
}
