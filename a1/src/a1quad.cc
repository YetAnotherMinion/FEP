#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>

#include <stdio.h>
#include <iostream>
  
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, false);
  //begin my code
  int x_elms = 4;
  int y_elms = 6;
  int total_elms = (x_elms+1) * (y_elms+1);
  double x_size = 8.0;
  double y_size = 8.0;
  //create an empty numbering
  apf::Numbering* numbers = apf::createNumbering(
      m,"my_numbers", m->getShape(), 1);
  int node_number = 0;
  //create an array to hold all our vertices
  apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
  for(int counter = 0; counter < total_elms; counter++) {
    vertices[counter] = m->createVert(0);
  }
  //create a pointer to pass in quad vertices
  apf::MeshEntity* quad_verts[4];

  //rather than test for that zero,zero element every time
  //just instantiate it manually once
  apf::Vector3* temp_vec = new apf::Vector3();
  double xyz[3] = {0,0,0};
  temp_vec->fromArray(xyz);
  int vert_index = 0; // use this to place verts in right spot
  m->setPoint(vertices[vert_index],0,*(temp_vec));
  //numbering each node in order of creation
  apf::number(numbers,vertices[vert_index],0,0, ++node_number);
  std::cout << node_number << ": " << *(temp_vec) << std::endl;
  
  for(int y_c = 1; y_c <= y_elms; y_c++) {
    //for the first column of x we make element to left
    xyz[0] = 0;
    xyz[1] = static_cast<double>(y_c) * y_size / 
      static_cast<double> (y_elms);
    temp_vec->fromArray(xyz); //load the values
    vert_index = (x_elms+1) * y_c;
    //print the location of each vertex
    m->setPoint(vertices[vert_index],0,*(temp_vec));
    apf::number(numbers,vertices[vert_index],0,0, ++node_number);
    std::cout << node_number << ": " << *(temp_vec) << std::endl;

    for(int x_c = 1; x_c <= x_elms; x_c++) {
      xyz[0] = static_cast<double>(x_c) * x_size /
	  static_cast<double>(x_elms);
      //for the first row of y we construct an extra element below
      if(y_c == 1) {
	xyz[1] = 0;
	temp_vec->fromArray(xyz);
	vert_index = x_c;
	m->setPoint(vertices[vert_index],0,*(temp_vec));
	apf::number(numbers,vertices[vert_index],0,0,++node_number);
	std::cout << node_number << ": " << *(temp_vec) << std::endl;
      } //now create the actual row
      xyz[1] = static_cast<double>(y_c) * y_size / 
	static_cast<double>(y_elms);
      temp_vec->fromArray(xyz);
      vert_index = y_c*(x_elms+1) + x_c;
      m->setPoint(vertices[vert_index],0,*(temp_vec));
      apf::number(numbers,vertices[vert_index],0,0,++node_number);
      //print the location of the vertex
      std::cout << node_number << ": " << *(temp_vec) << std::endl;
      //create the quad element
      quad_verts[0] = vertices[vert_index];
      quad_verts[1] = vertices[(vert_index - 1)];
      quad_verts[2] = vertices[(vert_index - x_elms - 2)];
      quad_verts[3] = vertices[(vert_index - x_elms - 1)];
      apf::buildElement(m, 0, apf::Mesh::QUAD, quad_verts);
    }
  }
  apf::deriveMdsModel(m);//this makes CAD model for classification
  //accept the changes
  m->acceptChanges();
  m->verify();
  
  //now look for getting all elements
  apf::Numbering* num_elms = numberOwnedDimension(m, "num_elms", 2);
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Downward down;
  int numVertices;
  apf::Vector3 point;
  int elm_num;
  while((e = m->iterate(it))) {
    /*get the vertices on this face element*/
    numVertices = m->getDownward(e, 0, down);
    elm_num = apf::getNumber(num_elms, e, 0, 0);
    std::cout << elm_num << ": ";
    for(int index = 0; index < numVertices; index++) {
      int vertex_id = apf::getNumber(numbers, down[index],0 ,0);
      std::cout << vertex_id << ' ';
      }
    std::cout << std::endl;
  }
  m->end(it);
  //provided: after this should not change
  apf::writeASCIIVtkFiles("outQuad", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
