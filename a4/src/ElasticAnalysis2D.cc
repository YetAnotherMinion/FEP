#include "ElasticAnalysis2D.h"
#include "StiffnessContributor2D.h"
#include "ForceContributor2D.h"
#include "MeshAdjReorder.h"

#include <fstream>

#define NUM_COMPONENTS 2

apf::Vector3 dummy(apf::Vector3 const& p){
	return apf::Vector3(10,0,0);
}


ElasticAnalysis2D::ElasticAnalysis2D(struct ElasticAnalysisInput & in)
{
	this->integration_order = in.integration_order;
	this->m = in.m;
	/*we create a 3 dimensional field but only will use 2 dimensions of that, (x & y)*/
	this->field = createField(this->m, "Field_1", apf::VECTOR, this->m->getShape());
	apf::zeroField(this->field);
	/*find the Lame parameters for plane strain*/
	double lambda, mu;
	lambda = (in.Nu * in.E) / ((1.0 + in.Nu) * (1.0 - 2.0 * in.Nu));
	mu = in.E / (2.0 + 2.0 * in.Nu);
	/*modify for plane stress*/
	const uint32_t PLANE_STRESS = 1;
	if(PLANE_STRESS) {
		lambda = (2.0 * lambda * mu) / ( lambda + 2.0 * mu);
	}
	/*initialize the D matrix for which ever case we are using*/
	this->D[0][0] = lambda + 2 * mu;
	this->D[0][1] = lambda;
	this->D[0][2] = 0.0;
	this->D[1][0] = lambda;
	this->D[1][1] = lambda + 2 * mu;
	this->D[1][2] = 0.0;
	this->D[2][0] = 0.0;
	this->D[2][1] = 0.0;
	this->D[2][2] = mu;
	/*set up element numberings we will use for assembly*/
	this->nodeNums = apf::createNumbering(this->m, "nodeNums", this->m->getShape(), NUM_COMPONENTS);
    this->faceNums = apf::createNumbering(this->m, "faceNums", apf::getConstant(this->m->getDimension()), 1);
	if(in.reorder == true){
		adjReorder(this->m, this->m->getShape(), NUM_COMPONENTS, this->nodeNums, this->faceNums);
		/*TEMPORARY: view the mesh*/
		apf::writeVtkFiles("elasticQuad", this->m);
	}
	/*compute the global degrees of freedom for the mesh*/
	std::size_t n_global_dofs = apf::countNodes(this->nodeNums) * NUM_COMPONENTS;
	/*initialize the linear system*/
	this->linsys = new AlgebraicSystem(n_global_dofs);
}

ElasticAnalysis2D::~ElasticAnalysis2D()
{
	delete this->linsys;
}

uint32_t ElasticAnalysis2D::setup()
{
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*iterate over the faces first*/
	it = this->m->begin(2);
	while((e = this->m->iterate(it))) {
		this->makeStiffnessContributor(e);
		this->makeForceContributor(e);
	}
	this->m->end(it);
	/*then pick up the edges*/
	it = this->m->begin(1);
	while((e = this->m->iterate(it))) {
		//this->makeForceContributor(e);
	}
	return 0;
}

uint32_t ElasticAnalysis2D::solve()
{
	std::ofstream out_file;
	out_file.open("K.txt", std::ios::out);
	out_file << this->linsys->K << std::endl;
	out_file.close();
	return 0;
}

uint32_t ElasticAnalysis2D::makeStiffnessContributor(apf::MeshEntity* e)
{
	int entity_type = this->m->getType(e);
	/*stiffness contributors are entirely local to this method, and are
	* assembled directly as they are generated. For meshes of same element
	* type the allocation will not change from element to element. For mixed
	* meshes, the interface will remain the same, but the allocation size will
	* be stuck on the largest element size. For example if the largest element
	* is a 4th order Lagrange element then the size of the local stiffness
	* contributor will grow as the entities are processed up to the size of
	* the largest element, and after that it will remain this size until the
	* end of the program. This is acceptable right now since we do
	* intend to do mixed meshes with different orders. Even if they are
	* the same order, the difference in size between different elements
	* will increase as the order of the element shape functions increase.
	* We will stick to low order (1 or 2) shape functions so this size
	* descrapancy is minimal. For comparison see table below for number of
	* lagrange nodes of Nth order for triangular and quadrilateral elements
	* order 	1 	2 	3	4
	* tri  		3 	6   10  15
	* quad 		4   9	16	25
	* As one can see. the difference between a triangular and quadrilateral
	* element is a difference of ((50 dofs)^2 - (30 dofs)^2) is allocating
	* 78% more storage than is actually required for processing the triangular
	* element. Clearly a strong canidate for causing unneccesary cache
	* misses.
	**/
	StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

	if(entity_type == apf::Mesh::QUAD || entity_type == apf::Mesh::TRIANGLE){
		apf::MeshElement* me = apf::createMeshElement(this->m, e);
		stiff.process(me);
		apf::destroyMeshElement(me);
		/*view the intermediate matrix*/
		uint32_t ndofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
		apf::NewArray< int > node_mapping(ndofs);
		apf::getElementNumbers(nodeNums, e, node_mapping);
		// std::cout  << "length of numbering: " << length << std::endl;
		// for(uint32_t ii = 0; ii < nnodes; ++ii){
		// 	std::cout << "Node " << ii << ": " << node_mapping[ii] << std::endl;
		// }
		this->linsys->assemble(stiff.ke, node_mapping, ndofs);

	} else {
		/*only accepts faces, so indicate improper input*/
		std::cout << "entity type: " << entity_type << std::endl;
		return 1;
	}
	return 0;
}

uint32_t ElasticAnalysis2D::makeForceContributor(apf::MeshEntity* e)
{
	int entity_type = this->m->getType(e);
	std::cout << "=====force contributor " << entity_type << " =======" << std::endl;

	apf::Vector3(*fnc_ptr)(apf::Vector3 const& p);
	fnc_ptr = &dummy;

	ForceContributor2D force(this->field, this->integration_order, fnc_ptr);

	apf::MeshElement* me = apf::createMeshElement(this->m, e);
	force.process(me);
	apf::destroyMeshElement(me);

	/*view the intermediate matrix*/
	uint32_t ndofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
	apf::NewArray< int > node_mapping(ndofs);
	apf::getElementNumbers(nodeNums, e, node_mapping);

	for(std::size_t ii = 0; ii < force.fe.size(); ++ii){
		std::cout << "Node " << node_mapping[ii] << ": " << force.fe[ii] << std::endl;
	}
	std::cout << "========================" << std::endl;
	// std::cout  << "length of numbering: " << length << std::endl;
	// for(uint32_t ii = 0; ii < nnodes; ++ii){
	// 	
	// }
	this->linsys->assemble(force.fe, node_mapping, ndofs);

	return 0;
}

uint32_t ElasticAnalysis2D::makeConstraint(apf::MeshEntity* e)
{
	return 0;
}

uint32_t ElasticAnalysis2D::recover()
{
	return 0;
}
