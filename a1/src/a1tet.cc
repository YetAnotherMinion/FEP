#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <stdio.h>
#include <iostream>

double computeTetVolume(const apf::Vector3 vec[]);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "tet-mesh-1.smb");
  apf::writeVtkFiles("beforeOutTet", m);
  //
  // insert mesh query code here
  //
  std::cout << "Part 1" << std::endl;
  apf::Numbering* region_nums = apf::numberOwnedDimension(m, "vol", 3);
  apf::Numbering* vert_nums = apf::numberOwnedDimension(m, "verts", 0);
  apf::Numbering* edge_nums = apf::numberOwnedDimension(m, "edges", 1);
  apf::MeshEntity* e;
  apf::Vector3 points[4];
  apf::MeshIterator* it = m->begin(3);
  while((e = m->iterate(it))){
    apf::Downward vertices;
    int num_vertices = m->getDownward(e, 0, vertices);
    int region_id = apf::getNumber(region_nums, e, 0, 0);
    std::cout << region_id << ": ";
    for(int i = 0; i < num_vertices; i++) {
      m->getPoint(vertices[i], 0 , points[i]);
      int vert_id = apf::getNumber(vert_nums,vertices[i],0 ,0);
      std::cout << vert_id << " ";
    }
    std::cout << std::endl;
    //compute the volume of each tet
    double volume = computeTetVolume((points));
    std::cout << "\tvolume: " << volume << std::endl;
  }
  m->end(it);
  //start #2
  std::cout << std::endl << "Part 2" << std::endl;
  apf::Adjacent up;
  it = m->begin(0);
  while((e = m->iterate(it))) {
    int vert_id = apf::getNumber(vert_nums, e, 0, 0);
    m->getAdjacent(e, 1, up);
    int num_edges = up.getSize();
    std::cout << vert_id << ": ";
    for(int i = 0; i < num_edges ; ++i) {
      int edge_id = apf::getNumber(edge_nums, up[i], 0 ,0);
      std::cout << edge_id << " ";
    }
    std::cout << std::endl;
  }
  m->end(it);
  std::cout << std::endl << "Part 6" << std::endl;
  //start #6, change the mesh to quadratic
  apf::changeMeshShape(m, apf::getLagrange(2));
  apf::Vector3 point;
  m->acceptChanges();
  m->verify();
  //now print out a list of each mid-edge node
  it = m->begin(1);//get the edges
  while((e = m->iterate(it))) {;
    m->getPoint(e, 0, point);
    std::cout << point << std::endl;
  }
  m->end(it);
  std::cout << std::endl;
  m->acceptChanges();
  m->verify();
  apf::writeVtkFiles("afterOutTet", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

double computeTetVolume(const apf::Vector3 vec[4]) {
  double array[4][4];
  for(int i = 0; i < 4; ++i) {
    array[i][0] = vec[i].x();
    array[i][1] = vec[i].y();
    array[i][2] = vec[i].z();
    array[i][3] = 1;
  }
  apf::Matrix<4,4> total(array);
  double result = apf::getDeterminant(total);
  result = result < 0 ? result*-1.0 : result;
  return result /6.0;
}
