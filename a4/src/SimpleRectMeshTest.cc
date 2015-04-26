#include <stdint.h>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>

#include "MeshBuilder.h"


class RectMeshTest : public testing::Test
{
protected:
	int five;
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
		/*use mesh builder to abstract mesh creations*/
		mesh_builder = new MeshBuilder();
		printf("virtual void setup is called\n");
		mesh = NULL;
		
		printf("%x\n", (long) mesh );
		mesh_builder->build2DRectQuadMesh(mesh, 40, 50, 10, 10, -10, -20);
		five = 6;
		printf("%x\n", (long) mesh );
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
			printf("Tear down called\n");
		}	
		//delete mesh_builder;
	}
	/*Helper method*/
	static int timesSeven(int n) {
		return n * 7;
	}
};

TEST_F(RectMeshTest, FirstTest) {
	printf("%d\n", five );
	EXPECT_EQ(7, timesSeven(1));
	apf::writeVtkFiles("outQuad", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("second", mesh);
}

TEST_F(RectMeshTest, SecondTest) {
	EXPECT_EQ(14, timesSeven(2));
}