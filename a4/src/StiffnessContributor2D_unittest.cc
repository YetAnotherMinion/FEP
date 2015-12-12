#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <cmath>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>

#include "StiffnessContributor2D.h"
#include "MeshBuilder.h"
#include "ElasticAnalysis2D.h"
#include "MeshAdjReorder.h" 

enum MeshTypes {
	LINEAR_QUAD,
	QUADRATIC_QUAD,
	LINEAR_TRI,
	QUADRATIC_TRI,
	SERENDIPITY_QUAD
};

struct test_parameters_wrapper {
	test_parameters_wrapper(MeshTypes mt, int intorder) 
		: mesh_type(mt), integration_order(intorder) {}
	MeshTypes mesh_type;
	int integration_order;
};

class StiffnessTest : public testing::Test,
	public ::testing::WithParamInterface<struct test_parameters_wrapper>
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;
	apf::Field* field;
	apf::Matrix< 3,3 > D;
	double E;
	double Nu;
	int integration_order;

	void changeMeshFromIndex(MeshTypes index) {
		if(NULL != this->mesh) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		int X_ELMS = 2;
		int Y_ELMS = 1;
		switch(index) {
			case LINEAR_QUAD:
				mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(1));
				break;
			case QUADRATIC_QUAD:
				this->mesh_builder->build2DRectQuadMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(2));
				break;
			case LINEAR_TRI:
				this->mesh_builder->build2DRectTriMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				break;
			case QUADRATIC_TRI:
				this->mesh_builder->build2DRectTriMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(2));
				break;
			case SERENDIPITY_QUAD:
				this->mesh_builder->build2DRectQuadMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getSerendipity());
			default:
				/*mesh should stay null*/
				break;
		}
		bool use_plane_stress = true;
		this->D = buildD(this->E, this->Nu, use_plane_stress);
		this->field = createField(this->mesh, "dummy", apf::VECTOR, this->mesh->getShape());
		apf::zeroField(this->field); /*must zero the field to force writing of data*/

	}

	virtual void SetUp() {
		mesh_builder = new MeshBuilder();
		mesh = NULL;
		this->E = 1e8;
		this->Nu = 0.35;
		this->integration_order = 4;
		struct test_parameters_wrapper tmp_foo(LINEAR_QUAD, 7);
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		delete mesh_builder;	
	}

};

/*C++ allows us to leave off the "struct" in a typename*/
INSTANTIATE_TEST_CASE_P(DifferentMeshOrders, StiffnessTest,
	::testing::Values(  test_parameters_wrapper(LINEAR_QUAD, 2),
					  	test_parameters_wrapper(LINEAR_QUAD,3),
					  	test_parameters_wrapper(LINEAR_QUAD,4),
					  	test_parameters_wrapper(QUADRATIC_QUAD,3),
					  	test_parameters_wrapper(QUADRATIC_QUAD,4),
					  	test_parameters_wrapper(LINEAR_TRI,2),
					  	test_parameters_wrapper(LINEAR_TRI,3),
					  	test_parameters_wrapper(LINEAR_TRI,4),
					  	test_parameters_wrapper(QUADRATIC_TRI,3),
					  	test_parameters_wrapper(QUADRATIC_TRI,4),
					  	test_parameters_wrapper(SERENDIPITY_QUAD,3),
					  	test_parameters_wrapper(SERENDIPITY_QUAD,4)));

TEST_P(StiffnessTest, StiffnessIsSymmetric) {
	/*use a linear quad*/
	struct test_parameters_wrapper tmp_args = GetParam();
	/*change the integration order away from the default of 4*/
	this->integration_order = tmp_args.integration_order;
	changeMeshFromIndex(tmp_args.mesh_type);

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = this->mesh->begin(2);
	while((e = this->mesh->iterate(it))) {
		/*compute each of the nodal submatrices by different method, and
		* do not explictly take advantage of symmetry*/
		StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

		apf::MeshElement* me = apf::createMeshElement(this->mesh, e);
		stiff.process(me);
		apf::destroyMeshElement(me);

		apf::Element* f_elm = apf::createElement(this->field, me);
		uint32_t nnodes = apf::countNodes(f_elm);
		apf::destroyElement(f_elm);

		uint32_t n_eqs = nnodes * this->mesh->getDimension();
		/*check that stiffness are symmetric*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = ii+1; jj < n_eqs; ++jj) {
				EXPECT_FLOAT_EQ(stiff.ke(ii,jj), stiff.ke(jj,ii));
			}
		}
	}
	this->mesh->end(it);

}

TEST_P(StiffnessTest, CheckStiffnessMatrix) {
	/*use a linear quad*/
	struct test_parameters_wrapper tmp_args = GetParam();
	/*change the integration order away from the default of 4*/
	this->integration_order = tmp_args.integration_order;
	changeMeshFromIndex(tmp_args.mesh_type);

	/*use a small local D, as any error in calculation betwee the
	* two methods will be linearly increased with magnitude of D*/
	this->D[0][0] = 1;
	this->D[0][1] = 1;
	this->D[1][0] = 1;
	this->D[1][1] = 1;
	this->D[2][2] = 1;


	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = this->mesh->begin(2);
	while((e = this->mesh->iterate(it))) {
		/*compute each of the nodal submatrices by different method, and
		* do not explictly take advantage of symmetry*/

		StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

		apf::MeshElement* me = apf::createMeshElement(this->mesh, e);
		stiff.process(me);

		apf::Element* f_elm = apf::createElement(this->field, me);
		uint32_t nnodes = apf::countNodes(f_elm);
		uint32_t n_eqs = nnodes * this->mesh->getDimension();

		apf::DynamicMatrix alt_ke;

		alt_ke.setSize(n_eqs,n_eqs);
		alt_ke.zero();

		int entity_type = this->mesh->getType(e);
		uint32_t nIntPoints = apf::countGaussPoints(entity_type, this->integration_order);
		/*compute the differential volume once, then construct a new strain tensor*/
		double jac_det = apf::getDV(me, apf::Vector3(0.0, 0.0, 0.0));
		apf::Matrix< 3,3 > local_D = this->D * jac_det;

		for(uint32_t ii = 0; ii < nIntPoints; ++ii) {
			apf::Vector3 p(0.0,0.0,0.0);
			apf::getGaussPoint(entity_type, this->integration_order, ii, p);
			double weight = apf::getIntWeight(me, this->integration_order, ii);

			apf::NewArray<apf::Vector3> gradShape;
			apf::getShapeGrads(f_elm, p, gradShape);

			for(uint32_t A = 0; A < nnodes; ++A) {
				for(uint32_t B = 0; B < nnodes; ++B) {
					/* nodal_stiffness_11
					* in reduced form we have N_A,x * D_11 * N_B,x
											+ N_A,y * D_33 * N_B,y*/
					double tmp = gradShape[A][0] * gradShape[B][0] * local_D[0][0];
					tmp += gradShape[A][1] * gradShape[B][1] *  local_D[2][2];
					tmp *= weight;
					alt_ke(2*A,2*B) += tmp;
					/* nodal_stiffness_12
					* we have N_A,x * D_12 * N_B,y
							+ N_A,y * D_33 * N_B,x*/
					tmp = gradShape[A][0] * gradShape[B][1] * local_D[0][1];
					tmp += gradShape[A][1] * gradShape[B][0] * local_D[2][2];
					tmp *= weight;
					alt_ke(2*A, 2*B + 1) += tmp;
					/* nodal_stiffness_21
					*we have N_A,y * D_12 * N_B,x
							+N_A,x * D_33 * N_B,y*/
					tmp = (gradShape[A][1] * gradShape[B][0]) * local_D[0][1];
					tmp += (gradShape[A][0] * gradShape[B][1])  * local_D[2][2];
					tmp *= weight;
					alt_ke(2*A + 1, 2*B) += tmp;
					/* nodal_stiffness_22
					*we have N_A,y * D_22 * N_B,y
							+N_A,x * D_33 * N_B,x*/
					tmp = gradShape[A][1] * gradShape[B][1] * local_D[1][1];
					tmp += gradShape[A][0] * gradShape[B][0] * local_D[2][2] ;
					tmp *= weight;
					alt_ke(2*A + 1, 2*B + 1) += tmp;
				}
			}
		}
		/*check that altenate calculation of stiffness is symmetric*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = ii+1; jj < n_eqs; ++jj) {
				EXPECT_FLOAT_EQ(alt_ke(ii,jj), alt_ke(jj,ii));
			}
		}

		/*now check that we computed the same stiffness matrix
		* by the two different methods*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = ii; jj < n_eqs; ++jj) {
				EXPECT_NEAR(alt_ke(ii,jj), stiff.ke(ii,jj), 3e-15);
			}
		}
		apf::destroyElement(f_elm);
		apf::destroyMeshElement(me);
	}
	this->mesh->end(it);
}

TEST_F(StiffnessTest, IdenticalLinearElement) {
	this->mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1,
					0.0, 0.0, 2.0, 1.0);
	apf::changeMeshShape(this->mesh, apf::getLagrange(1));
	this->field = createField(this->mesh, "dummy", apf::VECTOR, this->mesh->getShape());
	apf::zeroField(this->field); /*must zero the field to force writing of data*/
	EXPECT_TRUE(this->mesh != NULL);


	/*set up element numberings we will use for assembly*/
	uint32_t num_components = 2;
	apf::Numbering* nodeNums = apf::createNumbering(this->mesh, NODE_NUM_TAG_NAME, this->mesh->getShape(), num_components);
	apf::Numbering* faceNums = apf::createNumbering(this->mesh, FACE_NUM_TAG_NAME, apf::getConstant(this->mesh->getDimension()), 1);
	
	adjReorder(this->mesh, this->mesh->getShape(), num_components, nodeNums, faceNums);
	apf::writeASCIIVtkFiles("identical linear mesh", this->mesh);

	this->D[0][0] = 1;
	this->D[0][1] = 1;
	this->D[1][0] = 1;
	this->D[1][1] = 1;
	this->D[2][2] = 1;

	apf::MeshIterator* it;
	apf::MeshEntity* e1, *e2;
	it = this->mesh->begin(2);
	e1 = this->mesh->iterate(it);
	e2 = this->mesh->iterate(it);
	this->mesh->end(it);

	this->integration_order = 4;
	std::cout << this->D << std::endl;
	StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

	apf::DynamicMatrix ke1, ke2;
	apf::MeshElement* me;

	me = apf::createMeshElement(this->mesh, e1);
	int entity_type = this->mesh->getType(e1);
	uint32_t n_l_dofs = apf::countElementNodes(this->mesh->getShape(), entity_type) * num_components;
	apf::NewArray< int > node_mapping_1(n_l_dofs);
	int tmp_sz = apf::getElementNumbers(nodeNums, e1, node_mapping_1);
	assert(tmp_sz >= 0);
	/*we want to make sure that our guess for the vector was the right size
	* so we know how many degress of freedom our nodes have*/
	assert((uint32_t)tmp_sz == n_l_dofs);

	stiff.process(me);
	ke1 = stiff.ke;


	me = apf::createMeshElement(this->mesh, e2);
	entity_type = this->mesh->getType(e2);
	n_l_dofs = apf::countElementNodes(this->mesh->getShape(), entity_type) * num_components;
	apf::NewArray< int > node_mapping_2(n_l_dofs);
	tmp_sz = apf::getElementNumbers(nodeNums, e2, node_mapping_2);
	assert(tmp_sz >= 0);
	stiff.process(me);
	ke2 = stiff.ke;

	std::size_t nrows1, ncols1, nrows2, ncols2;
	nrows1 = ke1.getRows();
	ncols1 = ke1.getColumns();
	nrows2 = ke2.getRows();
	ncols2 = ke2.getColumns();

	assert(nrows1 == nrows2);
	assert(ncols1 == ncols2);

	std::cout << "ncols: " << ncols1 << " nrows: " << nrows1 << std::endl;

	for(std::size_t ii = 0; ii < nrows1; ++ii) {
		for(std::size_t jj = 0; jj < ncols1; ++jj) {
			EXPECT_FLOAT_EQ(ke1(ii, jj), ke2(ii, jj));
		}
	}
}


TEST_F(StiffnessTest, IdenticalQuadElement) {
	this->mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1,
					0.0, 0.0, 2.0, 1.0);
	apf::changeMeshShape(this->mesh, apf::getLagrange(2));
	this->field = createField(this->mesh, "dummy", apf::VECTOR, this->mesh->getShape());
	apf::zeroField(this->field); /*must zero the field to force writing of data*/
	EXPECT_TRUE(this->mesh != NULL);


	/*set up element numberings we will use for assembly*/
	uint32_t num_components = 2;
	apf::Numbering* nodeNums = apf::createNumbering(this->mesh, NODE_NUM_TAG_NAME, this->mesh->getShape(), num_components);
	apf::Numbering* faceNums = apf::createNumbering(this->mesh, FACE_NUM_TAG_NAME, apf::getConstant(this->mesh->getDimension()), 1);
	
	adjReorder(this->mesh, this->mesh->getShape(), num_components, nodeNums, faceNums);
	apf::writeASCIIVtkFiles("identical mesh", this->mesh);

	this->D[0][0] = 1;
	this->D[0][1] = 1;
	this->D[1][0] = 1;
	this->D[1][1] = 1;
	this->D[2][2] = 1;

	/*This is ordering of shape functions used below
	*      ^n
	*      |
	*   4--7--3
	*   |  |  |
	*   8  9--6--->e
	*   |     |
	*   1--5--2
	*
	*   indices are just above minus one
	*/

	/*manually set the degree of freedom mapping for the two elements
	* based on the node numberings determined by using the reordering
	* code and visually matching nodes with above orientation*/

	apf::MeshIterator* it;
	apf::MeshEntity* e1, *e2;
	it = this->mesh->begin(2);
	e1 = this->mesh->iterate(it);
	e2 = this->mesh->iterate(it);
	this->mesh->end(it);

	this->integration_order = 4;
	std::cout << this->D << std::endl;
	StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

	apf::DynamicMatrix ke1, ke2;
	apf::MeshElement* me;

	me = apf::createMeshElement(this->mesh, e1);
	int entity_type = this->mesh->getType(e1);
	uint32_t n_l_dofs = apf::countElementNodes(this->mesh->getShape(), entity_type) * num_components;
	apf::NewArray< int > node_mapping_1(n_l_dofs);
	int tmp_sz = apf::getElementNumbers(nodeNums, e1, node_mapping_1);
	assert(tmp_sz >= 0);
	/*we want to make sure that our guess for the vector was the right size
	* so we know how many degress of freedom our nodes have*/
	assert((uint32_t)tmp_sz == n_l_dofs);

	stiff.process(me);
	ke1 = stiff.ke;


	me = apf::createMeshElement(this->mesh, e2);
	entity_type = this->mesh->getType(e2);
	n_l_dofs = apf::countElementNodes(this->mesh->getShape(), entity_type) * num_components;
	apf::NewArray< int > node_mapping_2(n_l_dofs);
	tmp_sz = apf::getElementNumbers(nodeNums, e2, node_mapping_2);
	assert(tmp_sz >= 0);
	stiff.process(me);
	ke2 = stiff.ke;

	std::size_t nrows1, ncols1, nrows2, ncols2;
	nrows1 = ke1.getRows();
	ncols1 = ke1.getColumns();
	nrows2 = ke2.getRows();
	ncols2 = ke2.getColumns();

	assert(nrows1 == nrows2);
	assert(ncols1 == ncols2);

	std::cout << "ncols: " << ncols1 << " nrows: " << nrows1 << std::endl;

	for(std::size_t ii = 0; ii < nrows1; ++ii) {
		for(std::size_t jj = ii+1; jj < ncols1; ++jj) {
			float ratio = ke1(ii,jj) / ke2(ii, jj);
			if(fabs(ratio - 1.0) > 1e-6) {
				std::cout << "ratio: "<< ratio << " " << ii << "," << jj << std::endl;
			}
			
			// EXPECT_FLOAT_EQ(ke1(ii, jj), ke2(ii, jj));
		}
	}
}

class MyUserFunction : public apf::Function
{
public:
	MyUserFunction(apf::Field* f, apf::Mesh* m) : node_field(f), mesh(m) {}

	void eval(apf::MeshEntity* e, double *result) {
		if(NULL != node_field) {
			apf::MeshElement* me = apf::createMeshElement(this->mesh, e);
			/*find the location of the single node on entity*/
			apf::Vector3 tmp_vec;
			apf::getVector(this->node_field, e, 0, tmp_vec);

			apf::Vector3 x_vec;
			apf::mapLocalToGlobal(me, tmp_vec, x_vec);
			*result = (2*x_vec[0] + 5*x_vec[1]);
		} else {
			std::cout << "failed to set field" << std::endl;
			*result = -1.0;
		}
	}
protected:
	apf::Field* node_field;
	apf::Mesh* mesh;
};

TEST_F(StiffnessTest, ScratchPad) {
	apf::Mesh2* m = NULL;
	this->mesh_builder->buildBatmanElementMesh(m);

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = m->begin(2);
	e = m->iterate(it);
	m->end(it);

	apf::Field* node_f = apf::createField(m, "nodeField", apf::SCALAR, m->getShape());
	apf::zeroField(node_f);
	MyUserFunction test_fnc(node_f, m);
	apf::Field* test_f = apf::createUserField(m, "testingField", apf::SCALAR, m->getShape(), &test_fnc);

	apf::MeshElement* me = apf::createMeshElement(m, e);
	apf::Element* f_elm = apf::createElement(test_f, me);

	/*pick a random spot inside the element coordinates*/
	apf::Vector3 sample_points[8] =
	{apf::Vector3(-0.1241,    0.4889, 0.0),
	 apf::Vector3( 0.4897,    0.0347, 0.0),
	 apf::Vector3( 0.4090,    0.7269, 0.0),
	 apf::Vector3( 0.4172,   -0.3034, 0.0),
	 apf::Vector3( 0.6715,    0.2939, 0.0),
	 apf::Vector3(-0.2075,   -0.7873, 0.0),
	 apf::Vector3( 0.7172,    0.8884, 0.0),
	 apf::Vector3( 0.6302,   -0.1471, 0.0)};
	for(uint32_t ii = 0; ii < 8; ++ii) {
		apf::Vector3 tmp_grad(0.0, 0.0, 0.0);
		apf::getGrad(f_elm, sample_points[ii], tmp_grad);
		EXPECT_FLOAT_EQ(2.0, tmp_grad[0]);
		EXPECT_FLOAT_EQ(5.0, tmp_grad[1]);
	} 
	apf::writeASCIIVtkFiles("batman_elm", m);
}

TEST_F(StiffnessTest, QuadraticShapeFunctionCheck) {



}
