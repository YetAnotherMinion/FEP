#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <apfNumbering.h>

#include "MeshBuilder.h"


class RectMeshTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

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
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0, 0, 2, 1);
	//apf::writeVtkFiles("outQuad", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	//apf::changeMeshShape(mesh, apf::getLagrange(1));
	//apf::writeVtkFiles("secondQuad", mesh);
	apf::MeshEntity* e;
	apf::MeshIterator* it = mesh->begin(2);
	//apf::Field * master_f = createFieldOn(mesh, "foo", apf::VECTOR); 
	apf::Field* master_f = createFieldOn(mesh, "foo", apf::SCALAR);
	apf::FieldShape* fs = mesh->getShape();

	apf::Numbering* nodes_numbers = apf::numberOwnedNodes(mesh, "fatman", fs);

	apf::Numbering* all_node_nums = apf::createNumbering(mesh, "allNodes", mesh->getShape(), 1);
	

	std::cout << "outside" << std::endl;
	int dofs = apf::tempAdjReorder(all_node_nums, mesh, fs, 1);
	
	//apf::MeshTag* other_nodes = apf::reorder(mesh, "batman");

	apf::writeVtkFiles("secondQuad", mesh);
	
	while((e = mesh->iterate(it))) {
		//apf::MeshElement* mesh_elm = apf::getMeshElement()
		//apf::Element* element_ptr = apf::createElement(master_f, e );
		//apf::EntityShape::getValues
		apf::Vector3 param;
		double zeros[3] = {0, 0, 0};

		param.fromArray(zeros);
		int order = 4; //this is the order accuracy, not the polynomial degree
		apf::MeshElement* mesh_elm = apf::createMeshElement(mesh, e);
		int num_int_points = apf::countIntPoints(mesh_elm, order);
		apf::EntityShape* es = fs->getEntityShape(mesh->getType(e));
		uint32_t num_nodes = apf::countElementNodes(fs, mesh->getType(e));

		apf::NewArray< int > node_mapping(num_nodes);

		apf::getElementNumbers(nodes_numbers, e, node_mapping);

		for(int ii = 0; ii < num_nodes; ++ii) {
			std::cout << ii << " => " << node_mapping[ii] << std::endl;
		}

		for(int ii = 0; ii < num_int_points; ++ ii) {
			double running_sum = 0;

			apf::getIntPoint(mesh_elm, order, ii, param);
			std::cout << std::setprecision(16);
			std::cout << ii << ": " << param.x() << " " << param.y() << std::endl;

			apf::NewArray< apf::Vector3 > grads(num_nodes);
			es->getLocalGradients(mesh, e, param, grads);

			for(uint32_t jj = 0; jj < num_nodes; ++jj ) {
				std::cout << grads[jj].x() << " " << grads[jj].y() << " " << grads[jj].z() << std::endl;
				//apf::Vector3 xi(3,3,3);
				//std::cout << ii << ": " << xi.x() << ", " << xi.y() << ", " << xi.z() << std::endl;
				//fs->getNodeXi(apf::Mesh::QUAD, ii, xi);
				//std::cout << ii << ": " << xi.x() << ", " << xi.y() << ", " << xi.z() << std::endl;
			}
			std::cout << "========" << ii << "========" << std::endl;
		}
		std::cout << apf::countElementNodes(fs, mesh->getType(e)) << std::endl;
		//std::cout << master_f->countComponents() << std::endl;
	}

}

TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}
