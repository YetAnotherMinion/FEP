#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <PCU.h>

#include "MeshBuilder.h"

#include <unistd.h>

#define SECRET_BUILDER_NUMBERING "SecretBuilderNumbering"

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
  	PCU_Comm_Init();

	apf::Mesh2* mesh = NULL;
	MeshBuilder* mesh_builder = new MeshBuilder();
	mesh_builder->buildBatmanElementMesh(mesh);
	std::cout << "end of my code" << std::endl;
	usleep(1000000);

	apf::writeASCIIVtkFiles("batman", mesh);

	PCU_Comm_Free();
	MPI_Finalize();
	return 0;
}
