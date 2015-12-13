#include <stdint.h>
#include <assert.h>

#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "MeshBuilder.h"
#include "GeometryMappings.h"

EntityNumberer::EntityNumberer(apf::Mesh* m)
{
	this->mesh = m;
	this->tag = NULL;
	this->number_to_apply = NONE;
}

EntityNumberer::~EntityNumberer()
{

}

void EntityNumberer::call(apf::MeshEntity *e)
{
	if(this->tag != NULL) {
		this->mesh->setIntTag(e, this->tag, &this->number_to_apply);
	}
}

MeshBuilder::MeshBuilder() {
}
MeshBuilder::~MeshBuilder() {
}

void MeshBuilder::build2DRectQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
	uint32_t y_elms, double x0, double y0, double xf, double yf)
{
	/*there is nothing to be done for a zero size mesh so return*/
	if ((x_elms == 0) || (y_elms == 0)) {
		mesh = NULL;
		return;
	}
	gmi_register_null();
  	gmi_model* g = gmi_load(".null");
  	mesh = apf::makeEmptyMdsMesh(g, 2, false);

  	uint64_t total_elms = (x_elms + 1) * (y_elms + 1);
  	long double x_size = xf - x0;
  	long double y_size = yf - y0;

  	//create an empty numbering
	apf::Numbering* numbers = apf::createNumbering(mesh,SECRET_BUILDER_NUMBERING, 
												   mesh->getShape(), 1);

	apf::MeshTag* vertBCtag = mesh->createIntTag(VERT_BC_TAG_NAME, 1);
	apf::MeshTag* edgeBCtag = mesh->createIntTag(EDGE_BC_TAG_NAME, 1);
	apf::MeshTag* faceBCtag = mesh->createIntTag(FACE_BC_TAG_NAME, 1);

	EntityNumberer ent_number_functor(mesh);

	int node_number = 0;
	//create an array to hold all our vertices
	apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
	for(uint32_t counter = 0; counter < total_elms; counter++) {
		vertices[counter] = mesh->createVert(0);
	}
	//create a pointer to pass in quad vertices
	apf::MeshEntity* quad_verts[4];
	apf::MeshEntity* edge_verts[2];

	//rather than test for that zero,zero element every time
	//just instantiate it manually once
	apf::Vector3* temp_vec = new apf::Vector3();
	double xyz[3] = {x0, y0, 0};
	temp_vec->fromArray(xyz);
	int vert_index = 0; // use this to place verts in right spot
	mesh->setPoint(vertices[vert_index],0,*(temp_vec));
	//numbering each node in order of creation
	apf::number(numbers,vertices[vert_index],0,0, ++node_number);
	
	for(uint32_t y_c = 1; y_c <= y_elms; ++y_c) {
		//for the first column of x we make element to left
		xyz[0] = x0;
		xyz[1] = (static_cast<double>(y_c) * y_size / 
				  static_cast<double> (y_elms)) + y0;
		temp_vec->fromArray(xyz); //load the values
		vert_index = (x_elms + 1) * y_c;
		//print the location of each vertex
		mesh->setPoint(vertices[vert_index],0,*(temp_vec));
		apf::number(numbers,vertices[vert_index],0,0, ++node_number);

		for(uint32_t x_c = 1; x_c <= x_elms; x_c++) {
			xyz[0] = (static_cast<double>(x_c) * x_size / 
					  static_cast<double>(x_elms)) + x0;
			//for the first row of y we construct an extra element below
			if(y_c == 1) {
				xyz[1] = y0;
				temp_vec->fromArray(xyz);
				vert_index = x_c;
				mesh->setPoint(vertices[vert_index],0,*(temp_vec));
				apf::number(numbers,vertices[vert_index],0,0,++node_number);
			} //now create the actual row
			xyz[1] = (static_cast<double>(y_c) * y_size / 
					  static_cast<double>(y_elms)) + y0;
			temp_vec->fromArray(xyz);
			vert_index = y_c*(x_elms+1) + x_c;
			mesh->setPoint(vertices[vert_index],0,*(temp_vec));
			apf::number(numbers,vertices[vert_index],0,0,++node_number);
			/*create the quad element*/
			quad_verts[0] = vertices[vert_index];
			quad_verts[1] = vertices[(vert_index - 1)];
			quad_verts[2] = vertices[(vert_index - x_elms - 2)];
			quad_verts[3] = vertices[(vert_index - x_elms - 1)];
			/*with a one element mesh these conditions are not mutally exclusive*/
			if(x_c == x_elms) {
				/*the quad we are making is on the right edge*/
				if(y_c == y_elms) {
					/*this means we know where the x1y1 vertex is now*/
					int tmp_number_to_apply = VERT_X1Y1;
					// std::cout << "vert x1y1: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index)], vertBCtag, &tmp_number_to_apply);
				}
				if(1 == y_c) {
					int tmp_number_to_apply = VERT_X1Y0;
					// std::cout << "vert x1y0: " << apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index - x_elms -1)], vertBCtag, &tmp_number_to_apply);
				}
				// std::cout << "Right: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0) << ", " 
				// 	<< std::endl;

				ent_number_functor.number_to_apply = RIGHT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index)];
				edge_verts[1] = vertices[(vert_index - x_elms -1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);

			} 
			if( x_c == 1) {
				if(y_c == y_elms) {
					/*this means we know where the x1y1 vertex is now*/
					int tmp_number_to_apply = VERT_X0Y1;
					// std::cout << "vert x0y1: " << apf::getNumber(numbers, vertices[(vert_index-1)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index-1)], vertBCtag, &tmp_number_to_apply);
				}
				if(1 == y_c) {
					int tmp_number_to_apply = VERT_X0Y0;
					// std::cout << "vert x0y0: " << apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index- x_elms -2)], vertBCtag, &tmp_number_to_apply);
				}
				// std::cout << "Left: " << apf::getNumber(numbers, vertices[(vert_index-1)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = LEFT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index-1)];
				edge_verts[1] = vertices[(vert_index- x_elms -2)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			if(y_c == y_elms) {
				// std::cout << "Top: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index-1)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = TOP_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index)];
				edge_verts[1] = vertices[(vert_index-1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			if(1 == y_c) {
				// std::cout << "Bot: " << apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = BOT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index- x_elms -2)];
				edge_verts[1] = vertices[(vert_index - x_elms -1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			ent_number_functor.number_to_apply = ALL_FACES;
			ent_number_functor.tag = faceBCtag;
			apf::buildElement(mesh, NULL, apf::Mesh::QUAD, quad_verts, &ent_number_functor);
		}
	}
	apf::deriveMdsModel(mesh);/*this makes CAD model for classification*/
	/*accept the changes*/

	/*uncomment below to remove the debugging numbering used in the unittests
	* WARNING: do not destroy this if you are building the tests, as it will
	* fail a good number of those*/
	// apf::destroyNumbering(numbers);
	mesh->acceptChanges();
	mesh->verify();
	delete temp_vec;
	delete[] vertices;

}
void build2DJaggedQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf) 
{

	mesh = NULL;
}
void MeshBuilder::build2DRectTriMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
	uint32_t y_elms, double x0, double y0, double xf, double yf) 
{
	/*there is nothing to be done for a zero size mesh so return*/
	if ((x_elms == 0) || (y_elms == 0)) {
		mesh = NULL;
		return;
	}
	gmi_register_null();
  	gmi_model* g = gmi_load(".null");
  	mesh = apf::makeEmptyMdsMesh(g, 2, false);

  	uint64_t total_elms = (x_elms + 1) * (y_elms + 1);
  	long double x_size = xf - x0;
  	long double y_size = yf - y0;

  	//create an empty numbering
	apf::Numbering* numbers = apf::createNumbering(mesh,SECRET_BUILDER_NUMBERING, 
												   mesh->getShape(), 1);

	apf::MeshTag* vertBCtag = mesh->createIntTag(VERT_BC_TAG_NAME, 1);
	apf::MeshTag* edgeBCtag = mesh->createIntTag(EDGE_BC_TAG_NAME, 1);
	apf::MeshTag* faceBCtag = mesh->createIntTag(FACE_BC_TAG_NAME, 1);

	EntityNumberer ent_number_functor(mesh);
	int node_number = 0;
	//create an array to hold all our vertices
	apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
	for(uint32_t counter = 0; counter < total_elms; counter++) {
		vertices[counter] = mesh->createVert(0);
	}
	//create a pointer to pass in quad vertices
	apf::MeshEntity* tri_verts[3];
	apf::MeshEntity* edge_verts[2];

	//rather than test for that zero,zero element every time
	//just instantiate it manually once
	apf::Vector3* temp_vec = new apf::Vector3();
	double xyz[3] = {x0, y0, 0};
	temp_vec->fromArray(xyz);
	int vert_index = 0; // use this to place verts in right spot
	mesh->setPoint(vertices[vert_index],0,*(temp_vec));
	//numbering each node in order of creation
	apf::number(numbers,vertices[vert_index],0,0, ++node_number);
	
	for(uint32_t y_c = 1; y_c <= y_elms; ++y_c) {
		//for the first column of x we make element to left
		xyz[0] = x0;
		xyz[1] = (static_cast<double>(y_c) * y_size / 
				  static_cast<double> (y_elms)) + y0;
		temp_vec->fromArray(xyz); //load the values
		vert_index = (x_elms + 1) * y_c;
		//print the location of each vertex
		mesh->setPoint(vertices[vert_index],0,*(temp_vec));
		apf::number(numbers,vertices[vert_index],0,0, ++node_number);

		for(uint32_t x_c = 1; x_c <= x_elms; x_c++) {
			xyz[0] = (static_cast<double>(x_c) * x_size / 
					  static_cast<double>(x_elms)) + x0;
			//for the first row of y we construct an extra element below
			if(y_c == 1) {
				xyz[1] = y0;
				temp_vec->fromArray(xyz);
				vert_index = x_c;
				mesh->setPoint(vertices[vert_index],0,*(temp_vec));
				apf::number(numbers,vertices[vert_index],0,0,++node_number);
			} //now create the actual row
			xyz[1] = (static_cast<double>(y_c) * y_size / 
					  static_cast<double>(y_elms)) + y0;
			temp_vec->fromArray(xyz);
			vert_index = y_c*(x_elms+1) + x_c;
			mesh->setPoint(vertices[vert_index],0,*(temp_vec));
			apf::number(numbers,vertices[vert_index],0,0,++node_number);
			/*with a one element mesh these conditions are not mutally exclusive*/
			if(x_c == x_elms) {
				/*the quad we are making is on the right edge*/
				if(y_c == y_elms) {
					/*this means we know where the x1y1 vertex is now*/
					int tmp_number_to_apply = VERT_X1Y1;
					// std::cout << "vert x1y1: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index)], vertBCtag, &tmp_number_to_apply);
				}
				if(1 == y_c) {
					int tmp_number_to_apply = VERT_X1Y0;
					// std::cout << "vert x1y0: " << apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index - x_elms -1)], vertBCtag, &tmp_number_to_apply);
				}
				// std::cout << "Right: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0) << ", " 
				// 	<< std::endl;

				ent_number_functor.number_to_apply = RIGHT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index)];
				edge_verts[1] = vertices[(vert_index - x_elms -1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);

			} 
			if( x_c == 1) {
				if(y_c == y_elms) {
					/*this means we know where the x1y1 vertex is now*/
					int tmp_number_to_apply = VERT_X0Y1;
					// std::cout << "vert x0y1: " << apf::getNumber(numbers, vertices[(vert_index-1)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index-1)], vertBCtag, &tmp_number_to_apply);
				}
				if(1 == y_c) {
					int tmp_number_to_apply = VERT_X0Y0;
					// std::cout << "vert x0y0: " << apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0) << std::endl;

					mesh->setIntTag(vertices[(vert_index- x_elms -2)], vertBCtag, &tmp_number_to_apply);
				}
				// std::cout << "Left: " << apf::getNumber(numbers, vertices[(vert_index-1)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = LEFT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index-1)];
				edge_verts[1] = vertices[(vert_index- x_elms -2)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			if(y_c == y_elms) {
				// std::cout << "Top: " << apf::getNumber(numbers, vertices[(vert_index)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index-1)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = TOP_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index)];
				edge_verts[1] = vertices[(vert_index-1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			if(1 == y_c) {
				// std::cout << "Bot: " << apf::getNumber(numbers, vertices[(vert_index- x_elms -2)],0,0) << ", "
				// 	<< apf::getNumber(numbers, vertices[(vert_index - x_elms -1)],0,0)
				// 	<< std::endl;
				ent_number_functor.number_to_apply = BOT_EDGE;
				ent_number_functor.tag = edgeBCtag;
				edge_verts[0] = vertices[(vert_index- x_elms -2)];
				edge_verts[1] = vertices[(vert_index - x_elms -1)];
				apf::buildElement(mesh, NULL, apf::Mesh::EDGE, edge_verts, &ent_number_functor);
			}
			/*only need to set this up once since both triangles get same value*/
			ent_number_functor.number_to_apply = ALL_FACES;
			ent_number_functor.tag = faceBCtag;
			//creat two triangular elements from quad
			tri_verts[0] = vertices[vert_index];
			tri_verts[1] = vertices[(vert_index - 1)];
			tri_verts[2] = vertices[(vert_index - x_elms - 2)];
			apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, tri_verts, &ent_number_functor);
			tri_verts[0] = vertices[(vert_index - x_elms - 2)];
			tri_verts[1] = vertices[(vert_index - x_elms - 1)];
			tri_verts[2] = vertices[vert_index];
			apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, tri_verts, &ent_number_functor);
		}
	}
	apf::deriveMdsModel(mesh);//this makes CAD model for classification
	//accept the changes
	mesh->acceptChanges();
	mesh->verify();
	delete temp_vec;
	delete[] vertices;
}

void MeshBuilder::buildBatmanElementMesh(apf::Mesh2* & mesh)
{
		/*Batman element*/
	gmi_register_null();
	gmi_model* g = gmi_load(".null");
	mesh = apf::makeEmptyMdsMesh(g, 2, false);

	apf::Numbering* numbers = apf::createNumbering(mesh,SECRET_BUILDER_NUMBERING, 
												   mesh->getShape(), 1);
	apf::MeshEntity** vertices = new apf::MeshEntity*[4];
	for(uint32_t counter = 0; counter < 4; counter++) {
		vertices[counter] = mesh->createVert(0);
		apf::number(numbers,vertices[counter],0,0,counter);
	}


	/*this element is flipped upside down around x axis*/
	apf::Vector3 tmp_vec(1.35,2.36,0.0);
	mesh->setPoint(vertices[0], 0, tmp_vec);
	tmp_vec.x() = 2.57; tmp_vec.y() = 2.70;
	mesh->setPoint(vertices[1], 0, tmp_vec);
	tmp_vec.x() = 2.51; tmp_vec.y() = 1.75;
	mesh->setPoint(vertices[2], 0, tmp_vec);
	tmp_vec.x() = 1.64; tmp_vec.y() = 1.28;
	mesh->setPoint(vertices[3], 0, tmp_vec);

	apf::buildElement(mesh, NULL, apf::Mesh::QUAD, vertices);

	apf::deriveMdsModel(mesh);
	mesh->acceptChanges();
	mesh->verify();

	apf::changeMeshShape(mesh, apf::getLagrange(2));
	apf::MeshIterator *it;
	apf::MeshEntity* e;
	it = mesh->begin(2);
	e = mesh->iterate(it);

	apf::Downward down;
	uint32_t sz = mesh->getDownward(e, 1, down);
	assert(4 == sz);
	tmp_vec.x() = 2.09; tmp_vec.y() = 2.36;
	mesh->setPoint(down[0], 0, tmp_vec);
	tmp_vec.x() = 2.73; tmp_vec.y() = 2.25;
	mesh->setPoint(down[1], 0, tmp_vec);
	tmp_vec.x() = 2.01; tmp_vec.y() = 1.72;
	mesh->setPoint(down[2], 0, tmp_vec);
	tmp_vec.x() = 1.15; tmp_vec.y() = 1.77;
	mesh->setPoint(down[3], 0, tmp_vec);
	/*set the center node*/
	tmp_vec.x() = 1.90; tmp_vec.y() = 2.05;
	mesh->setPoint(e, 0, tmp_vec);
	/*clean up the iterator*/
	mesh->end(it);
	delete[] vertices;
}
