#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfNumbering.h>

#include "MeshBuilder.h"

#include <unistd.h>

#define SECRET_BUILDER_NUMBERING "SecretBuilderNumbering"

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
  	PCU_Comm_Init();

	apf::Mesh2* mesh = NULL;
	MeshBuilder* mesh_builder = new MeshBuilder();
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 1, 1);

	

	std::cout << "end of my code" << std::endl;
	usleep(1000000);

	apf::writeVtkFiles("tri_test", mesh);

	PCU_Comm_Free();
	MPI_Finalize();

	return 0;
}