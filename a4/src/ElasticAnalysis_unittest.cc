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

#include "TestingUtilityFunctions.h"

#define YOUNGS_MODULUS  1e8
#define POISSONS_RATIO 0.35
#define PRINT_STRESS_AND_STRAIN 0

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

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());

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

class SampleProblems : public ::testing::Test,
	public ::testing::WithParamInterface<MeshTypes> 
{
protected:
	apf::Mesh2* mesh;
	std::string suffix;

	virtual void SetUp() {
		mesh = NULL;
		suffix = "";
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
	}
};

INSTANTIATE_TEST_CASE_P(GeometryTags, SampleProblems,
	::testing::Values(MeshTypes::LINEAR_QUAD, MeshTypes::QUAD_QUAD,
					  MeshTypes::LINEAR_TRI, MeshTypes::QUAD_TRI,
					  MeshTypes::SERENDIPITY_QUAD));


TEST_P(SampleProblems, NonZeroDirchlet) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 10;
	uint32_t Y_ELMS = 10;
	ASSERT_TRUE(NULL == this->mesh);
	this->mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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

	/*Non zero dirchlet boundary condition*/
	cnstr_ptr = &smallDisplacementNegativeX_2D;
	geo_map->addDircheletMapping(RIGHT_EDGE, cnstr_ptr);
	
	struct ElasticAnalysisInput input = {
			this->mesh,
			geo_map,
			integration_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
		// MatView(tmp.linsys->K, PETSC_VIEWER_STDOUT_WORLD);

	//VecView(tmp.linsys->F, PETSC_VIEWER_STDOUT_WORLD);

	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

	std::cout << "************************************" << std::endl;
	std::cout << "*\t" << this->suffix << std::endl;
	std::cout << "************************************" << std::endl;
	std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	std::cout << std::endl;
	std::string out_name = "dirchlet";
	out_name += this->suffix;
	
	// apf::writeASCIIVtkFiles(out_name.c_str(), this->mesh);
	apf::writeVtkFiles(out_name.c_str(), this->mesh);
	delete geo_map;
}


TEST_P(SampleProblems, Gravity) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 10;
	uint32_t Y_ELMS = 10;
	ASSERT_TRUE(NULL == this->mesh);
	this->mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

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
	uint32_t X_ELMS = 10;
	uint32_t Y_ELMS = 10;
	this->mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

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

TEST_P(SampleProblems, ConstantTraction) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 2;
	uint32_t Y_ELMS = 2;
	this->mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

	std::cout << "************************************" << std::endl;
	std::cout << "*\t\t" << this->suffix << std::endl;
	std::cout << "************************************" << std::endl;
	std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	std::cout << std::endl;
	std::string bar = "constant_load";
	bar += this->suffix;
	apf::writeASCIIVtkFiles(bar.c_str(), this->mesh);

	delete geo_map;
}

TEST_P(SampleProblems, QuadraticTraction) {
	MeshTypes index = GetParam();
	uint32_t X_ELMS = 10;
	uint32_t Y_ELMS = 10;
	this->mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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
	geo_map->addDircheletMapping(BOT_EDGE, cnstr_ptr);

	apf::Vector3 (*traction_ptr)(apf::Vector3 const &);
	traction_ptr = &QuadraticLoad_X;
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

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

	std::cout << "************************************" << std::endl;
	std::cout << "*\t\t" << this->suffix << std::endl;
	std::cout << "************************************" << std::endl;
	std::cout << "strain energy: " << std::setprecision(20) << tmp.strain_energy;
	std::cout << std::endl;
	std::string bar = "quadratic_load";
	bar += this->suffix;
	apf::writeASCIIVtkFiles(bar.c_str(), this->mesh);

	delete geo_map;
}

TEST_P(SampleProblems, StrainEnergyLinearBodyForce) {
	MeshTypes index = GetParam();

	uint32_t X_ELMS = 10;
	uint32_t Y_ELMS = 10;
	apf::Mesh2* mesh = getMeshFromIndex(index, X_ELMS, Y_ELMS, this->suffix);
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

#if PRINT_STRESS_AND_STRAIN
	for(auto kv : tmp.stress) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "sigma_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "sigma_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "sigma_xy = " << kv.second[2] << std::endl;
	}
	std::cout << std::endl;
	for(auto kv : tmp.strain) {
		std::cout << "x = " << kv.first << std::endl << "\t";
		std::cout << "e_xx = " << kv.second[0] << std::endl << "\t";
		std::cout << "e_yy = " << kv.second[1] << std::endl << "\t";
		std::cout << "e_xy = " << kv.second[2] << std::endl;
	}
#endif

	/*preserve a record of some of these meshes*/
	std::string out_name("linear_body_force_2x2");
	out_name += this->suffix;
	apf::writeASCIIVtkFiles(out_name.c_str(), mesh);


	delete tmp;
	mesh->destroyNative();
	apf::destroyMesh(mesh);
	delete geo_map;
}