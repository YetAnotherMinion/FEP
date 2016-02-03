#include "ElasticAnalysis2D.h"
#include "StiffnessContributor2D.h"
#include "ForceContributor2D.h"
#include "MeshAdjReorder.h" 
#include "MeshBuilder.h"/*gives us the entity tags used find constraints
						* and boundary conditions */
#include "GeometryMappings.h" /*gives us the special enity tags enum*/
#include "RecoverAtIntegrationPoints.h"

#include <fstream>
#include <assert.h>
#include <stdexcept> /*provides runtime error*/

#include <PCU.h>
#include <apf.h>

#define NUM_COMPONENTS 2

apf::Vector3 dummy(apf::Vector3 const& p){
	return apf::Vector3(10,0,0);
}

apf::Matrix< 3,3 > buildD(double E, double Nu, bool use_plane_stress) {
	apf::Matrix< 3,3 > D;
	/*find the Lame parameters for plane strain*/
	double lambda, mu;
	lambda = (Nu * E) / ((1.0 + Nu) * (1.0 - 2.0 * Nu));
	mu = E / (2.0 + 2.0 * Nu);
	/*modify for plane stress*/
	if(use_plane_stress) {
		lambda = (2.0 * lambda * mu) / ( lambda + 2.0 * mu);
	}
	/*initialize the D matrix for which ever case we are using*/
	D[0][0] = lambda + 2 * mu;
	D[0][1] = lambda;
	D[0][2] = 0.0;
	D[1][0] = lambda;
	D[1][1] = lambda + 2 * mu;
	D[1][2] = 0.0;
	D[2][0] = 0.0;
	D[2][1] = 0.0;
	D[2][2] = mu;
	return D;
}

ElasticAnalysis2D::ElasticAnalysis2D(struct ElasticAnalysisInput & in) :
	strain_energy(0.0),
	geometry_map(in.geo_map),
	integration_order(in.integration_order),
	m(in.m)
{
	/*we create a 3 dimensional field but only will use 2 dimensions of that, (x & y)*/
	this->field = createField(this->m, "dummyField", apf::VECTOR, this->m->getShape());
	apf::zeroField(this->field);

	/*thisfield will store the resulting displacments*/
	this->disp_field = apf::createFieldOn(this->m, DISPLACEMENT_FIELD_NAME, apf::VECTOR);
	apf::zeroField(this->disp_field);

	bool use_plane_stress = true;
	this->D = buildD(in.E, in.Nu, use_plane_stress);

	/*set up element numberings we will use for assembly*/
	this->nodeNums = apf::createNumbering(this->m, NODE_NUM_TAG_NAME, this->m->getShape(), NUM_COMPONENTS);
	this->faceNums = apf::createNumbering(this->m, FACE_NUM_TAG_NAME, apf::getConstant(this->m->getDimension()), 1);
	if(in.reorder == true){
		adjReorder(this->m, this->m->getShape(), NUM_COMPONENTS, this->nodeNums, this->faceNums);
	}
	/*compute the global degrees of freedom for the mesh*/
	std::size_t n_global_dofs = apf::countNodes(this->nodeNums) * NUM_COMPONENTS;
	/*initialize the linear system*/
	this->linsys = new AlgebraicSystem(n_global_dofs);
	this->nGlobalDOFs = n_global_dofs;
}

ElasticAnalysis2D::~ElasticAnalysis2D()
{
	apf::destroyNumbering(this->nodeNums);
	apf::destroyNumbering(this->faceNums);
	apf::destroyField(this->disp_field);
	delete this->linsys;
}

uint32_t ElasticAnalysis2D::setup()
{
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*find all of the boundary conditions fixed,
	* chech all edges classified on model edge and
	* all vertices classified on model vertex*/
	apf::MeshTag* face_tag = this->m->findTag(FACE_BC_TAG_NAME);
	apf::MeshTag* edge_tag = this->m->findTag(EDGE_BC_TAG_NAME);
	apf::MeshTag* vert_tag = this->m->findTag(VERT_BC_TAG_NAME);
	/*this should exist to use geomery based definition*/
	assert(NULL != face_tag);
	assert(NULL != edge_tag);
	assert(NULL != vert_tag);

	it = this->m->begin(1);
	while((e = this->m->iterate(it))) {
		this->makeConstraint(e);
	}
	this->m->end(it);
	it = this->m->begin(0);
	while((e = this->m->iterate(it))) {
		this->makeConstraint(e);
	}
	this->m->end(it);
	/*now that all fixed dofs have been accounted for,
	* allow assembly of force and stiffness contributors*/
	this->linsys->beginAssembly();

	/*iterate over the faces*/
	it = this->m->begin(2);
	while((e = this->m->iterate(it))) {
		this->makeStiffnessContributor(e);
		/*chech for body forces*/
		this->makeForceContributor(e);	
	}
	this->m->end(it);
	/*then pick up the edges*/
	it = this->m->begin(1);
	while((e = this->m->iterate(it))) {
		this->makeForceContributor(e);
	}
	this->m->end(it);

	/*now synchronize the linear system in preparation for solving*/
	this->linsys->synchronize();
	return 0;
}

uint32_t ElasticAnalysis2D::solve()
{
	// std::fstream myfile;
	// myfile.open("K.txt");
	// myfile << this->linsys->K;
	// myfile.close();

	double t0 = PCU_Time();
	this->linsys->solve();
	double t1 = PCU_Time();
	std::cout << "System solved in: " << t1- t0 << " seconds" <<std::endl;
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
	* element. Clearly a strong candidate for causing unneccesary cache
	* misses.
	**/
	if(entity_type == apf::Mesh::QUAD || entity_type == apf::Mesh::TRIANGLE){

		StiffnessContributor2D stiff(this->field, this->D, this->integration_order);
		apf::MeshElement* me = apf::createMeshElement(this->m, e);
		stiff.process(me);
		apf::destroyMeshElement(me);
		/*view the intermediate matrix*/
		uint32_t n_l_dofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
		apf::NewArray< int > node_mapping(n_l_dofs);
		int tmp_sz = apf::getElementNumbers(this->nodeNums, e, node_mapping);
		assert(tmp_sz >= 0);
		/*we want to make sure that our guess for the vector was the right size
		* so we know how many degress of freedom our nodes have*/
		assert((uint32_t)tmp_sz == n_l_dofs);
		this->linsys->assemble(stiff.ke, node_mapping, tmp_sz);
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
	apf::Vector3(*fnc_ptr)(apf::Vector3 const &);
	fnc_ptr = NULL;
	/*there will be only one tag per element we consider*/
	apf::MeshTag *mesh_tag;
	int tag_data;

	if(entity_type == apf::Mesh::QUAD || entity_type == apf::Mesh::TRIANGLE) {
		/*lookup a specific tagname assumed to correspond to tractions
		* via the GeometryMappings object*/
		mesh_tag = this->m->findTag(FACE_BC_TAG_NAME);
	} else if(entity_type == apf::Mesh::EDGE) {
		/*lookup a specific tagname assumed to correspond to tractions
		* via the GeometryMappings object*/
		mesh_tag = this->m->findTag(EDGE_BC_TAG_NAME);
	} else if(entity_type == apf::Mesh::VERTEX) {
		/*lookup a specific tagname assumed to correspond to tractions
		* via the GeometryMappings object*/
		mesh_tag = this->m->findTag(VERT_BC_TAG_NAME);
	} else {
		/*bad entity*/
		return 1;
	}
	if(false == this->m->hasTag(e, mesh_tag)){
		/*nothing to be done if there is no traction*/
		return 0;
	}
	/*if execution reaches this point after all of the returns above,
	* then we know there is a tag on this entity*/
	this->m->getIntTag(e, mesh_tag, &tag_data);
	if(1 == this->geometry_map->neumann_map.count(tag_data)) {
		/*if the mapping does exist, then retrieve the function
		* pointer that represents the loadind over the domain*/
		fnc_ptr = this->geometry_map->neumann_map[tag_data];
		ForceContributor2D force(this->field, this->integration_order, fnc_ptr);
		/*below three lines are the simplest expression of a very complicated
		* numerical integration routine*/
		apf::MeshElement* me = apf::createMeshElement(this->m, e);
		force.process(me);
		apf::destroyMeshElement(me);
		/*compute the mapping of the element dofs into global dofs*/
		uint32_t n_l_dofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
		apf::NewArray< int > node_mapping(n_l_dofs);
		int tmp_sz = apf::getElementNumbers(this->nodeNums, e, node_mapping);
		assert(tmp_sz >= 0);
		assert((uint32_t)tmp_sz == n_l_dofs);
		this->linsys->assemble(force.fe, node_mapping, tmp_sz);
	}
	return 0;
}

uint32_t ElasticAnalysis2D::makeConstraint(apf::MeshEntity* e)
{
	int entity_type = this->m->getType(e);
	apf::MeshTag *mesh_tag;
	int tag_data;

	if(entity_type == apf::Mesh::EDGE) {
		/*lookup a specific tagname assumed to correspond to tractions
		* via the GeometryMappings object*/
		mesh_tag = this->m->findTag(EDGE_BC_TAG_NAME);
	} else if(entity_type == apf::Mesh::VERTEX) {
		/*lookup a specific tagname assumed to correspond to tractions
		* via the GeometryMappings object*/
		mesh_tag = this->m->findTag(VERT_BC_TAG_NAME);
	} else {
		/*bad entity*/
		return 1;
	}
	if(false == this->m->hasTag(e, mesh_tag)){
		/*nothing to be done if there is no traction*/
		return 0;
	}
	/*if execution reaches this point after all of the returns above,
	* then we know there is a tag on this entity*/
	this->m->getIntTag(e, mesh_tag, &tag_data);
	if(1 == this->geometry_map->dirchelet_map.count(tag_data))
	{
		// std::cout << "found key" << std::endl;
		/*retrieve the specifc boundary contion from the store and evaluate it*/
		void(*fnc_ptr)(apf::MeshEntity*,
					   apf::Mesh*,
					   apf::Numbering*,
					   std::vector< uint64_t >&,
					   std::vector < double > & );
		fnc_ptr = this->geometry_map->dirchelet_map[tag_data];
		std::vector< uint64_t > fixed_mapping(0);
		std::vector< double > displacement(0);
		fnc_ptr(e, this->m, this->nodeNums, fixed_mapping, displacement);

		this->linsys->addBoundaryConstraint(displacement, fixed_mapping);
		
	} else {
		// std::cout << "key not found" << std::endl;
	}
	return 0;
}

uint32_t ElasticAnalysis2D::recover()
{
	this->displacement.clear();
	this->strain.clear();
	this->stress.clear();

	this->linsys->extractDisplacement(this->displacement);
	assert(this->displacement.size() == this->nGlobalDOFs);
	/*set the values of the solution field by iterating over the
	* faces. We will overwrite values several times by this method
	* but we consider the loss of efficiency acceptable to avoid
	* some tedious bookkeeping. ENHANCEMENT Can always add this 
	* later */
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*reuse the same vector since the third componenet never changes
	* for this 2D problem*/
	double tmp_arry[3] = {0.0, 0.0, 0.0};
	assert(NUM_COMPONENTS <= 3);
	apf::FieldShape* fs = apf::getShape(this->disp_field);

	for(std::size_t dim = 0; dim <= 2; ++dim) {
		if(fs->hasNodesIn(dim)) {
			it = this->m->begin(dim);
			while((e = this->m->iterate(it))) {
				int entity_type = this->m->getType(e);
				// std::cout << "entity type is " << entity_type << std::endl;
				int n_elm_nodes = fs->countNodesOn(entity_type);
				// std::cout << "Dim: "<<  dim << " There are N=" << n_elm_nodes << " nodes" << std::endl;
				assert (n_elm_nodes >= 0);
				for(std::size_t ii = 0; ii < static_cast<std::size_t>(n_elm_nodes); ++ii) {
					/*get each component*/
					for(std::size_t jj = 0; jj < NUM_COMPONENTS; ++jj) {
						if(apf::isNumbered(this->nodeNums, e, ii, jj)) {
							uint64_t gdof_indx = apf::getNumber(this->nodeNums, e, ii, jj);
							assert(0 <= gdof_indx);
							assert(this->nGlobalDOFs > gdof_indx);
							tmp_arry[jj] = this->displacement[gdof_indx];
						} else {
							throw std::runtime_error("Node not numbered");
						}
						
					}
					apf::Vector3 tmp_vec(tmp_arry);
					apf::setVector(this->disp_field, e, ii, tmp_vec);
				}
			}
			this->m->end(it);
		}
	}
	/*to find the strain energy we compute first rhs = Ku, then
	* multiply by transpose(u) * rhs which is just the dot product.
	* Then we multiply by 1/2. First we will compute the*/


	/*Stiffness contributors for 2D mesh are only the faces*/
	it = this->m->begin(2);
	while((e = this->m->iterate(it))){
		/*fill in integrator*/
		/*extract the local displacements for this element*/
		int entity_type = this->m->getType(e);
		int n_elm_nodes = apf::countElementNodes(this->m->getShape(), entity_type);
		uint32_t n_l_dofs = n_elm_nodes * NUM_COMPONENTS;
		apf::NewArray< int > node_mapping(n_l_dofs);
		int tmp_sz = apf::getElementNumbers(this->nodeNums, e, node_mapping);
		assert(tmp_sz >= 0);
		assert((uint32_t)tmp_sz == n_l_dofs);

		/*squash global dofs into per element dofs*/
		std::vector< apf::Vector<2> > nodal_displacements(n_elm_nodes);
		apf::DynamicVector elm_displacements(n_l_dofs);
		/*copy the local tid*/
		/*use the convention [x0, y0, x1, y1, x2, ...] for the order
		* of the local dofs in the node_mapping vector*/
		for(uint32_t ii = 0; ii < static_cast<uint32_t>(n_elm_nodes); ++ii) {
			double tmp_components[NUM_COMPONENTS];
			for(uint32_t jj = 0; jj < NUM_COMPONENTS; ++jj) {
				uint64_t local_element_index = ii * NUM_COMPONENTS + jj;
				int global_dof = node_mapping[local_element_index];
				/*copy over element displacements according to the local DOFs*/
				elm_displacements[local_element_index] = this->displacement[global_dof];
				/*and copy the same displacement into a apf::Vector for each
				* node*/
				tmp_components[jj] = this->displacement[global_dof];
			}
			nodal_displacements[ii] = apf::Vector<NUM_COMPONENTS>(tmp_components);
		}

		RecoverAtIntegrationPoints<3, 2> recovery_helper(
			this->field,
			this->D,
			nodal_displacements,
			this->integration_order,
			&this->strain,
			&this->stress);

		StiffnessContributor2D elm_stiffness(this->field,
											 this->D,
											 this->integration_order);
		apf::MeshElement* me = apf::createMeshElement(this->m, e);
		recovery_helper.process(me);
		elm_stiffness.process(me);
		/*compute part of the strain energy Kd */
		apf::DynamicVector tmp_rhs(n_l_dofs);
		assert(n_l_dofs == elm_stiffness.getNumberOfDOFs());
		

		for(uint32_t ii = 0; ii < n_l_dofs; ++ii) {
			apf::DynamicVector ke_row;
			elm_stiffness.ke.getRow(static_cast<std::size_t>(ii), ke_row);
			tmp_rhs[ii] = ke_row * elm_displacements;
		}
		double element_energy = tmp_rhs * elm_displacements;
		this->strain_energy += element_energy;

		apf::destroyMeshElement(me);
	}
	this->m->end(it);
	/*Compute*/
	this->strain_energy *= 0.5;

	std::cout << "==============Recovered=============" << std::endl;

	return 0;
}
