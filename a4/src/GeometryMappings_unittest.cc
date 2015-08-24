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

#include "MeshAdjReorder.h"
#include "ElasticAnalysis2D.h"
#include "MeshBuilder.h"
#include "GeometryMappings.h"
#include "AlgebraicSystem.h" /*provides SOLVER_ABSOLUTE_TOLERANCE constant*/

#define CONSTANT_DISPLACEMENT 3.5
void ConstantDisplacementX_2D(
	apf::MeshElement *me,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> & dofs,
	std::vector<double> & disp)
{
	int entity_type = me->getMesh()->getType(me->getEntity());
	uint32_t nnodes = apf::countElementNodes(me->getMesh()->getShape(), entity_type);
	/*in 2D we assume there are only 2 dofs per node*/
	apf::NewArray< int > node_mapping(nnodes*2);
	uint32_t tmp_sz = apf::getElementNumbers(nodeNums, me->getEntity(), node_mapping);
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
			disp[curs_indx] = CONSTANT_DISPLACEMENT;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}


class GeometryMappingTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	apf::Numbering* nodeNums;
	apf::Numbering* faceNums;

	virtual void SetUp()
	{
		this->mesh = NULL;
		MeshBuilder builder;
		builder.build2DRectQuadMesh(this->mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);

		this->nodeNums = apf::createNumbering(this->mesh, NODE_NUM_TAG_NAME, this->mesh->getShape(), 2);
		this->faceNums = apf::createNumbering(this->mesh, FACE_NUM_TAG_NAME, apf::getConstant(this->mesh->getDimension()), 1);
		adjReorder(this->mesh, this->mesh->getShape(), 2, this->nodeNums, this->faceNums);
	}
	virtual void TearDown()
	{
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
	}
};

TEST_F(GeometryMappingTest, NoneConstraint) {
	/*we do need to generate a mesh for this to work properly*/
	GeometryMappings geo_map;
	EXPECT_EQ(0, geo_map.dirchelet_map.count(LEFT_EDGE));
	
	void (*fnc_ptr)(apf::MeshElement*, apf::Numbering*, std::vector<uint64_t> &, std::vector<double> &);
	fnc_ptr = &noConstraint;

	geo_map.addDircheletMapping(LEFT_EDGE, fnc_ptr);
	EXPECT_EQ(1, geo_map.dirchelet_map.count(LEFT_EDGE));
}
