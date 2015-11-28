#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <set>
#include <algorithm>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <apfMatrix.h>

#include "MeshBuilder.h"
#include "MeshAdjReorder.h"
#include "GeometryMappings.h"


#define LINEAR_QUAD 0
#define QUAD_QUAD 1
#define LINEAR_TRI 2
#define QUAD_TRI 3
#define SERENDIPITY_QUAD 4

class RectMeshTest : public ::testing::Test,
	public ::testing::WithParamInterface<int> 
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	apf::Mesh2* getMeshFromIndex(int index, int X_ELMS, int Y_ELMS) {
		apf::Mesh2* tmp = NULL;
		switch(index) {
			case LINEAR_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				break;
			case QUAD_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getLagrange(2));
				break;
			case LINEAR_TRI:
				this->mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				break;
			case QUAD_TRI:
				this->mesh_builder->build2DRectTriMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getLagrange(2));
				break;
			case SERENDIPITY_QUAD:
				this->mesh_builder->build2DRectQuadMesh(tmp, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(tmp, apf::getSerendipity());
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
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
		delete mesh_builder;
	}
};

TEST_F(RectMeshTest, Rectangle) {
	int X_ELMS = 2;
	int Y_ELMS = 1;
	EXPECT_NO_THROW(mesh_builder->build2DRectQuadMesh(mesh, X_ELMS, Y_ELMS, 0.0, 0.0, 2.0, 2.0));
	EXPECT_TRUE(mesh != NULL);
}

TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	//apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getLagrange(2));
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}


INSTANTIATE_TEST_CASE_P(GeometryTags, RectMeshTest,
	::testing::Values(LINEAR_QUAD, QUAD_QUAD, LINEAR_TRI,
						QUAD_TRI, SERENDIPITY_QUAD));

TEST_P(RectMeshTest, CheckBoundaryEdgeTags) {
	int index = GetParam();
	uint32_t X_ELMS = 4;
	uint32_t Y_ELMS = 3;
	apf::Mesh2* specific_mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
	ASSERT_TRUE(specific_mesh != NULL);
	/*we do a bad thing and use undocumented numbering
	* that the mesh_builder creates to verify that numberings
	* are in the right place*/
	apf::MeshTag* edge_tag = specific_mesh->findTag(EDGE_BC_TAG_NAME);
	EXPECT_TRUE(NULL != edge_tag);
	apf::Numbering* secret_numbering = specific_mesh->findNumbering(SECRET_BUILDER_NUMBERING);
	EXPECT_TRUE(secret_numbering != NULL);
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*generate sets of the node numberings for the linear specific_mesh*/
	std::set<int32_t> left_edge;
	std::set<int32_t> right_edge;
	std::set<int32_t> top_edge;
	std::set<int32_t> bot_edge;
	uint32_t counter_left_edge = 0;
	uint32_t counter_right_edge = 0;
	uint32_t counter_top_edge = 0;
	uint32_t counter_bot_edge = 0;

	/*since there must always be at least one element in Y, and
	* there is special numbering behaviour for the bottom left 
	* element */
	left_edge.insert(1);
	left_edge.insert(2);
	right_edge.insert((2*X_ELMS + 1));
	for(std::size_t ii = 2*(X_ELMS+1); ii <= ((X_ELMS+1) * (Y_ELMS+1)); ii+=(X_ELMS+1)) {
		right_edge.insert(ii);
	}
	for(std::size_t ii = 1; ii <= (2*X_ELMS +1); ii+=2) {
		bot_edge.insert(ii);
	}
	for(std::size_t ii = (2* X_ELMS + 3); ii <= ((X_ELMS+1)* Y_ELMS +1); ii+=(X_ELMS+1)) {
		left_edge.insert(ii);
	}
	for(std::size_t ii = ((X_ELMS+1)* Y_ELMS +1); ii <= ((X_ELMS+1) * (Y_ELMS+1)); ++ii) {
		top_edge.insert(ii);
	}

	it = specific_mesh->begin(1);
	while((e = specific_mesh->iterate(it))) {
		apf::Downward down;
		int data;
		if(specific_mesh->hasTag(e, edge_tag)) {
			/*find the vertices bounding that edge*/
			int32_t n_verts = specific_mesh->getDownward(e, 0, down);
			EXPECT_TRUE(n_verts == 2);
			std::set<int> edge_verts;
			edge_verts.clear();
			assert(n_verts >= 0);
			for(std::size_t ii = 0; ii < static_cast<uint32_t>(n_verts); ++ii) {
				/*create a set of the edge vertices, this
				* really should only be two elements*/
				edge_verts.insert(apf::getNumber(secret_numbering, down[ii], 0,0));
			}
			std::set<int> intersect;
			intersect.clear();
			/*test left side*/
			set_intersection(left_edge.begin(), left_edge.end(), 
				edge_verts.begin(), edge_verts.end(), 
				std::inserter(intersect, intersect.begin()));
			if(intersect.size() == static_cast<uint32_t>(n_verts)) {
				EXPECT_TRUE(specific_mesh->hasTag(e, edge_tag));
				if(specific_mesh->hasTag(e, edge_tag)) {
					specific_mesh->getIntTag(e, edge_tag, &data);
					EXPECT_EQ(LEFT_EDGE, data);
					++counter_left_edge;
				}
			}
			/*check right edge*/
			intersect.clear();
			set_intersection(right_edge.begin(), right_edge.end(), 
				edge_verts.begin(), edge_verts.end(), 
				std::inserter(intersect, intersect.begin()));
			if(intersect.size() == static_cast<uint32_t>(n_verts)) {
				EXPECT_TRUE(specific_mesh->hasTag(e, edge_tag));
				if(specific_mesh->hasTag(e, edge_tag)) {
					specific_mesh->getIntTag(e, edge_tag, &data);
					EXPECT_EQ(RIGHT_EDGE, data);
					++counter_right_edge;
				}
			}
			/*check bot edge*/
			intersect.clear();
			set_intersection(bot_edge.begin(), bot_edge.end(), 
				edge_verts.begin(), edge_verts.end(), 
				std::inserter(intersect, intersect.begin()));
			if(intersect.size() == static_cast<uint32_t>(n_verts)) {
				EXPECT_TRUE(specific_mesh->hasTag(e, edge_tag));
				if(specific_mesh->hasTag(e, edge_tag)) {
					specific_mesh->getIntTag(e, edge_tag, &data);
					EXPECT_EQ(BOT_EDGE, data);
					++counter_bot_edge;
				}
			}
			/*check top edge*/
			intersect.clear();
			set_intersection(top_edge.begin(), top_edge.end(), 
				edge_verts.begin(), edge_verts.end(), 
				std::inserter(intersect, intersect.begin()));
			if(intersect.size() == static_cast<uint32_t>(n_verts)) {
				EXPECT_TRUE(specific_mesh->hasTag(e, edge_tag));
				if(specific_mesh->hasTag(e, edge_tag)) {
					specific_mesh->getIntTag(e, edge_tag, &data);
					EXPECT_EQ(TOP_EDGE, data);
					++counter_top_edge;
				}
			}
		}
	}
	specific_mesh->end(it);
	/*make sure that we see appropriate number of edges*/
	EXPECT_EQ(X_ELMS, counter_top_edge);
	EXPECT_EQ(X_ELMS, counter_bot_edge);
	EXPECT_EQ(Y_ELMS, counter_right_edge);
	EXPECT_EQ(Y_ELMS, counter_left_edge);

	/*clean up the generated mesh*/
	specific_mesh->destroyNative();
	apf::destroyMesh(specific_mesh);	
}

TEST_P(RectMeshTest, CheckBoundaryVertexTags) {
	int index = GetParam();
	int X_ELMS = 3;
	int Y_ELMS = 3;

	apf::Mesh2* specific_mesh = this->getMeshFromIndex(index, X_ELMS, Y_ELMS);
	ASSERT_TRUE(specific_mesh != NULL);
	/*we do a bad thing and use undocumented numbering
	* that the mesh_builder creates to verify that numberings
	* are in the right place*/
	apf::MeshTag* vert_tag = specific_mesh->findTag(VERT_BC_TAG_NAME);
	EXPECT_TRUE(NULL != vert_tag);
	apf::Numbering* secret_numbering = specific_mesh->findNumbering(SECRET_BUILDER_NUMBERING);
	EXPECT_TRUE(NULL != secret_numbering);
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = specific_mesh->begin(0);
	while((e = specific_mesh->iterate(it))) {
		int data;
		if(specific_mesh->hasTag(e, vert_tag)) {
			specific_mesh->getIntTag(e, vert_tag, &data);
			int node_label = apf::getNumber(secret_numbering, e, 0,0);
			/*we make sure that each corner of the rectangle is identifiable*/
			if( ((X_ELMS+1)*(Y_ELMS+1))== node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X1Y1) == data);
			} else if(1 == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X0Y0) == data);
			} else if((X_ELMS * (Y_ELMS + 1) +1) == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X0Y1) == data);
			} else if((X_ELMS * 2 + 1) == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X1Y0) == data);
			}
		}
	}
	specific_mesh->end(it);
	/*clean up the generated mesh*/
	specific_mesh->destroyNative();
	apf::destroyMesh(specific_mesh);	
}
