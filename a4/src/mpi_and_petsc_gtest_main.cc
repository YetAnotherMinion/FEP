#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <petscsys.h>

#include <gtest/gtest.h>

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	/*MPI init comes after google testing, since google test will remove
	* its specific arguments from argv*/
	MPI_Init(&argc,&argv);
  	PCU_Comm_Init();
  	PetscInitialize(0,0,0,0);

	int rc = RUN_ALL_TESTS();
	if(rc != 0){
		printf("Google Test returned: %d", rc);
	}
	
	PetscFinalize();
	PCU_Comm_Free();
  	MPI_Finalize();
	return 0;
}
