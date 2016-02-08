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

class ConvergenceProblems : public ::testing::Test,
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

INSTANTIATE_TEST_CASE_P(GeometryTags, ConvergenceProblems,
	::testing::Values(MeshTypes::LINEAR_QUAD, MeshTypes::QUAD_QUAD,
					  MeshTypes::LINEAR_TRI, MeshTypes::QUAD_TRI,
					  MeshTypes::SERENDIPITY_QUAD));

#define UPPER_H_LIMIT 80

TEST_P(ConvergenceProblems, StrainEnergyConstantBodyForce) {
	MeshTypes index = GetParam();

	double energies[UPPER_H_LIMIT];
	int problem_dofs[UPPER_H_LIMIT];

	for(uint32_t h = 1; h < UPPER_H_LIMIT; ++h) {
		uint32_t X_ELMS = h;
		uint32_t Y_ELMS = h+1;
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
		// if( 10 == h ) {
		// 	std::string out_name("constant_body_force_2x2");
		// 	out_name += this->suffix;
		// 	apf::writeASCIIVtkFiles(out_name.c_str(), mesh);
		// }

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

TEST_P(ConvergenceProblems, StrainEnergyLinearBodyForce) {
	MeshTypes index = GetParam();

	double energies[UPPER_H_LIMIT];
	int problem_dofs[UPPER_H_LIMIT];

	for(uint32_t h = 1; h < UPPER_H_LIMIT; ++h) {
		uint32_t X_ELMS = h;
		uint32_t Y_ELMS = h+1;
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

		problem_dofs[h] = tmp->linsys->freeDOFs;
		energies[h] = tmp->strain_energy;
		/*preserve a record of some of these meshes*/
		// if( 10 == h ) {
		// 	std::string out_name("linear_body_force_2x2");
		// 	out_name += this->suffix;
		// 	apf::writeASCIIVtkFiles(out_name.c_str(), mesh);
		// }


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

TEST_P(ConvergenceProblems, StrainEnergyQuadraticTraction) {
	MeshTypes index = GetParam();

	double energies[UPPER_H_LIMIT];
	int problem_dofs[UPPER_H_LIMIT];

	for(uint32_t h = 1; h < UPPER_H_LIMIT; ++h) {
		uint32_t X_ELMS = h;
		uint32_t Y_ELMS = h+1;
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
		traction_ptr = &QuadraticLoad_X;
		geo_map->addNeumannMapping(RIGHT_EDGE, traction_ptr);

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
		// if( 10 == h ) {
		// 	std::string out_name("linear_body_force_2x2");
		// 	out_name += this->suffix;
		// 	apf::writeASCIIVtkFiles(out_name.c_str(), mesh);
		// }

		delete tmp;
		mesh->destroyNative();
		apf::destroyMesh(mesh);
		delete geo_map;
	}

	std::cout << "quadratic_traction" << this->suffix << " = [" << std::endl; 
	for(auto ii = 1; ii < UPPER_H_LIMIT; ++ii) {
		std::cout << ii << ", " << problem_dofs[ii] << ", ";
		std::cout << std::setprecision(20);
		std::cout << energies[ii] <<  std::setprecision(5) << std::endl;
	}
	std::cout << "];" << std::endl;
}