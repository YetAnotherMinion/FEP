#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <assert.h>
#define INT_TAG_LENGTH 3

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "parallelMesh.smb");
  //Part 1 and 2
  //iterate through each vertex and compute adjacency for each
  apf::Numbering* vert_nums = apf::numberOwnedDimension(m, "verts", 0);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;
  apf::MeshEntity* e;
  apf::Adjacent up;
  apf::MeshTag* local_adj_tag = m->createIntTag("all_adj", 3);
  //iterate over each mesh vertex and get all adjacencies
  int current_partition = PCU_Comm_Self();
  int efr_adj[INT_TAG_LENGTH];
  //use tag size to record  with number of bytes in our tag
  int TRANSFER_TAG_SIZE = sizeof(int) * INT_TAG_LENGTH;
  PCU_Comm_Begin();
  while((vert = m->iterate(it))) {
    for(int dim = 1; dim < 4; ++dim) {
      m->getAdjacent(vert, dim, up);
      int owned_dim_count = 0;
      for(unsigned int ii = 0; ii < up.getSize(); ++ii) {
  //only count adjacencies on the local part
  //this avoid double counting
  if(m->isOwned(up[ii])) {
    ++owned_dim_count;
  }
      }
      //offset by one
      int index = dim - 1;
      efr_adj[index] = owned_dim_count;
    }
    //now add the data to the tag
    m->setIntTag(vert, local_adj_tag, efr_adj);
    //now send the taged data to the other partitions
    //only send the data that other partitions want
    if( ! m->isOwned(vert)) {
      apf::Copies remotes;
      m->getRemotes(vert, remotes);
      for(apf::Copies::iterator cit = remotes.begin(); cit != remotes.end(); ++cit) {
        //send the vertex first
        int target_process = cit->first;
        if(target_process != m->getOwner(vert)){
          continue;
        }
        apf::MeshEntity* over_there_vert = cit->second;
        int rc = PCU_Comm_Pack(target_process, &over_there_vert, sizeof(apf::MeshEntity*));
        rc = PCU_Comm_Pack(target_process, efr_adj, TRANSFER_TAG_SIZE );
      }
    }
  }
  //we are done with local adjacencies
  PCU_Comm_Send();
  
  while(PCU_Comm_Receive()) {
    apf::MeshEntity* target_vert;
    PCU_Comm_Unpack(&target_vert, sizeof(apf::MeshEntity*));
    //parts should send adjacency to correct part number
    assert(m->isOwned(target_vert));
    int foreign_adj[INT_TAG_LENGTH];
    PCU_Comm_Unpack(foreign_adj ,TRANSFER_TAG_SIZE );
    //get the current adj tag and add
    int temp_tag[INT_TAG_LENGTH];
    m->getIntTag(target_vert, local_adj_tag, temp_tag);    
    //sum the current adjacencies and the local adjacencies
    for(int jj = 0; jj < INT_TAG_LENGTH; ++jj) {
      temp_tag[jj] += foreign_adj[jj];
    } 
    m->setIntTag(target_vert, local_adj_tag, temp_tag);
  }
  //now report the adjacencies for each vertex in the local part
  it = m->begin(0);
  while((vert = m->iterate(it))) {
    if(m->isOwned(vert)) {
      int temp_tag[INT_TAG_LENGTH];
      m->getIntTag(vert, local_adj_tag, temp_tag);
      printf("L: %d Vert num: %d, E: %d F: %d R: %d\n", current_partition, \
       apf::getNumber(vert_nums, vert, 0, 0), temp_tag[0], temp_tag[1], temp_tag[2]);
    }
  }
  //Part 3
  apf::writeVtkFiles("pre_migr", m);
  //migration plan should exist for each process
  apf::Migration* plan = new apf::Migration(m);
  //get the mesh regions that are on boundary
  it = m->begin(2);
  apf::Adjacent up_r;
  //count the number of faces on the partition model face
  int partition_face_cnt = 0;
  while((e = m->iterate(it))) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    if( ! m->isOwned(e)) {
      m->getAdjacent(e, 3, up_r); 
      assert(up_r.getSize() == 1);
      ++partition_face_cnt;
      if(PCU_Comm_Self() == 0) { 
        apf::MeshEntity* region = up_r[0];
        if(! plan->has(region)) {
           plan->send(region, m->getOwner(e));
        }
      }
    }
  }
  printf("Partition %d reports %d partition faces before migration\n", \
   PCU_Comm_Self() , partition_face_cnt);
  //now send all of the migration plan
  m->migrate(plan);
  //now compute the number of shared faces
  it = m->begin(2);
  //count the number of faces on the partition model face
  int aftr_partition_face_cnt = 0;
  while((e = m->iterate(it))) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    if( ! m->isOwned(e)) {
      m->getAdjacent(e, 3, up_r); 
      assert(up_r.getSize() == 1);
      ++aftr_partition_face_cnt;
    }
  }
  printf("Partition %d reports %d partition faces after migration\n", \
   PCU_Comm_Self() , aftr_partition_face_cnt);
  apf::writeVtkFiles("post_migr", m); 
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
