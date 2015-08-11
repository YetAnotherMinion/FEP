#ifndef ELASTIC_ANALYSIS
#define ELASTIC_ANALYSIS 1

#include <apf.h>
#include <apfMesh.h>

#include "FEAnalysis.h"

class ElasticAnalysis2D : FEAnalysis
{
public:
	ElasticAnalysis2D(apf::Mesh* m);
	virtual ~ElasticAnalysis2D();

	virtual uint32_t setup();
	virtual uint32_t solve();
	virtual uint32_t makeStiffnessContributor(apf::MeshEntity* e);
	virtual uint32_t makeForceContributor(apf::MeshEntity* e);
	virtual uint32_t makeConstraint(apf::MeshEntity* e);
	virtual uint32_t recover();

private:
	/* data */
	apf::Mesh* m;
};

#endif