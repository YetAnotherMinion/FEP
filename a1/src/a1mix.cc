#include <apf.h>
#include <gmi_mesh.h>
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
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "mixed-mesh-1.smb");
  apf::Numbering* face_nums = numberOwnedDimension(m, "faces", 2);
  apf::MeshEntity* e;
  apf::ModelEntity* me;
  //brute force find the upwards classification
  apf::MeshIterator* it = m->begin(2);
  #define MODEL_FACE 81
  std::cout << "\nFaces classified on model face " << MODEL_FACE << std::endl;
  while((e = m->iterate(it))) {
    me = m->toModel(e);
    int dimension = m->getModelType(me);
    int tag = m->getModelTag(me);
    if(tag == MODEL_FACE) {
      int face_id = apf::getNumber(face_nums, e, 0, 0);
      std::cout << face_id << " " ;
    }
  }
  std::cout << std::endl << std::endl;

  //find the type of region for each entity
  apf::Numbering* region_nums = numberOwnedDimension(m ,"regions", 3);;  
  it = m->begin(3);
  while((e = m->iterate(it))) {
    int region_id = apf::getNumber(region_nums, e, 0, 0);
    std::cout << region_id << ": ";
    int region_type = m->getType(e);
    switch(region_type){
    case(apf::Mesh::TET):
      std::cout << "TET";
      break;
    case(apf::Mesh::HEX):
      std::cout << "HEX";
      break;
    case(apf::Mesh::PRISM):
      std::cout << "PRISM";
      break;
    case(apf::Mesh::PYRAMID):
      std::cout << "PYRAMID";
      break;
    default:
      std::cout << "not a region";
    };
    std::cout << std::endl;
  }
  m->end(it);

  apf::writeASCIIVtkFiles("outMixed", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
