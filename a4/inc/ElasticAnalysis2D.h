#ifndef ELASTIC_ANALYSIS
#define ELASTIC_ANALYSIS 1

#include <apf.h>
#include <apfMesh.h>
#include <apfMatrix.h>
#include <apfNumbering.h>

#include "FEAnalysis.h"
#include "AlgebraicSystem.h"
#include "GeometryMappings.h"

#define NODE_NUM_TAG_NAME "nodeNums"
#define FACE_NUM_TAG_NAME "faceNums"

#define DISPLACEMENT_FIELD_NAME "DisplacementField"
#define STRAIN_FIELD_NAME "StrainField"
#define STRESS_FIELD_NAME "StressField"


struct ElasticAnalysisInput {
	apf::Mesh* m;
	GeometryMappings* geo_map;
	uint32_t integration_order;
	double E;
	double Nu;
	bool reorder;
};

class ElasticAnalysis2D : FEAnalysis
{
public:
	ElasticAnalysis2D(struct ElasticAnalysisInput & in);
	virtual ~ElasticAnalysis2D();

	virtual uint32_t setup();
	virtual uint32_t solve();
	virtual uint32_t makeStiffnessContributor(apf::MeshEntity* e);
	virtual uint32_t makeForceContributor(apf::MeshEntity* e);
	virtual uint32_t makeConstraint(apf::MeshEntity* e);
	virtual uint32_t recover();
	
	std::vector<double> displacement;
	std::vector<double> strain;
	std::vector<double> stress;
	double strain_energy;

	AlgebraicSystem* linsys;
private:
	GeometryMappings* geometry_map;
	uint32_t integration_order;
	apf::Mesh* m;
	apf::Field* field;
	apf::Field* disp_field;
	apf::Matrix< 3,3> D;
	apf::Numbering* nodeNums;
	apf::Numbering* faceNums;
	std::size_t nGlobalDOFs;

};

/*expose this internal method for testing purposes*/
apf::Matrix< 3,3 > buildD(double E, double Nu, bool use_plane_stress);

#endif
