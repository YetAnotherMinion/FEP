#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <stdexcept> /*provides runtime error*/

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

		apf::Element* f_elm = apf::createElement(this->field, me);
		uint32_t nnodes = apf::countNodes(f_elm);
		apf::destroyElement(f_elm);
		apf::destroyMeshElement(me);

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

/*used to configure the location of the mesh element*/
struct point_and_step {
	point_and_step(double x, double y, double offx, double offy,
				 double xsz, double ysz, uint32_t shpordr) : x0(x), y0(y), offset_x(offx),
				 offset_y(offy), x_sz(xsz), y_sz(ysz) 
		{ 
			std::cout << "shape order: " << shpordr << std::endl;
			if(shpordr != 1 && shpordr != 2) {
				throw std::runtime_error("Only shape function orders 1 and 2 are supported");
			} else {
				this->shape_order = shpordr;
			}
		};
	double x0;
	double y0;
	/*the offset of the second mesh created with a sigle element*/
	double offset_x;
	double offset_y;
	/*the sizes of both elements need to be the same*/
	double x_sz;
	double y_sz;
	uint32_t shape_order;
};

class LocalStiffnessMatrixTest : public testing::Test,
	public ::testing::WithParamInterface<struct point_and_step>
{
protected:
	apf::Mesh2* m1, *m2;
	MeshBuilder* mesh_builder;
	apf::Field* field1, *field2;
	apf::Numbering* nodeNums;
	apf::Matrix< 3,3 > D;
	double E;
	double Nu;
	uint32_t num_components;
	int integration_order;

	void createTwoDifferentMeshes(struct point_and_step ps) {
		/*forcibly clean up any left over meshes*/
		if(NULL != this->m1) {
			m1->destroyNative();
			apf::destroyMesh(m1);
		}
		if(NULL != this->m2) {
			m2->destroyNative();
			apf::destroyMesh(m2);
		}

		int X_ELMS, Y_ELMS;
		X_ELMS = 1;
		Y_ELMS = 1;

		this->mesh_builder->build2DRectQuadMesh(m1, X_ELMS, Y_ELMS,
						ps.x0,
						ps.y0,
						ps.x0 + ps.x_sz,
						ps.y0 + ps.y_sz);
		apf::changeMeshShape(m1, apf::getLagrange(ps.shape_order));
		
		this->field1 = createField(m1, "dummy1", apf::VECTOR, m1->getShape());
		apf::zeroField(this->field1); /*must zero the field to force writing of data*/
		EXPECT_TRUE(m1 != NULL);
		/*Build second single mesh element far away*/
		this->mesh_builder->build2DRectQuadMesh(m2, X_ELMS, Y_ELMS,
						(ps.x0 + ps.offset_x),
						(ps.y0 + ps.offset_y),
						(ps.x0 + ps.offset_x + ps.x_sz),
						(ps.y0 + ps.offset_y + ps.y_sz));

		apf::changeMeshShape(m2, apf::getLagrange(ps.shape_order));
		field2 = createField(m2, "dummy2", apf::VECTOR, m2->getShape());
		apf::zeroField(field2); /*must zero the field to force writing of data*/
		EXPECT_TRUE(m2 != NULL);

		bool use_plane_stress = true;
		this->D = buildD(this->E, this->Nu, use_plane_stress);
	}

	virtual void SetUp() {
		/*2D mesh only has two componenets that are labeled for each node*/
		this->num_components = 2;
		mesh_builder = new MeshBuilder();
		m1 = NULL;
		m2 = NULL;
		this->E = 1e8;
		this->Nu = 0.35;
		this->integration_order = 4;
	}

	virtual void TearDown() {
		if(m1 != NULL) {
			m1->destroyNative();
			apf::destroyMesh(m1);
		}
		if(m2 != NULL) {
			m2->destroyNative();
			apf::destroyMesh(m2);
		}
		delete mesh_builder;	
	}

};
/*only test points close by (~+-10 away) from each other because 
* the numerical calculation of gradients is sensitive to the location in
* global coordinate space, so two identical elements that are millions
* of units apart will have an absolute error rougly millions of times larger*/
INSTANTIATE_TEST_CASE_P(IdenticalElement, LocalStiffnessMatrixTest,
	::testing::Values(  point_and_step(10000000.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1), /*points far way*/
						point_and_step(10000000.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2),
						point_and_step(-10000000.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1), /*check negative quadrant*/
						point_and_step(-10000000.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2)));

TEST_P(LocalStiffnessMatrixTest, DifferentLocations) {
	/*use a linear quad*/
	struct point_and_step tmp_args = GetParam();
	createTwoDifferentMeshes(tmp_args);

	/*prevent scaling of error between different stiffness matricies by linear amount*/
	this->D[0][0] = 1;
	this->D[0][1] = 1;
	this->D[1][0] = 1;
	this->D[1][1] = 1;
	this->D[2][2] = 1;
	/*set up element numberings we will use for assembly*/
	uint32_t num_components = 2;

	apf::MeshIterator* it;
	apf::MeshEntity* e1, *e2;
	it = m1->begin(2);
	e1 = m1->iterate(it);
	m1->end(it);
	it = m2->begin(2);
	e2 = m2->iterate(it);
	m2->end(it);
	
	this->integration_order = 4;
	StiffnessContributor2D stiff1(field1, this->D, this->integration_order);

	apf::DynamicMatrix ke1, ke2, ke3;
	apf::MeshElement* me;

	me = apf::createMeshElement(m1, e1);
	stiff1.process(me);
	apf::destroyMeshElement(me);
	me = NULL;
	ke1 = stiff1.ke;

	StiffnessContributor2D stiff2(field2, this->D, this->integration_order);

	me = apf::createMeshElement(m2, e2);
	stiff2.process(me);
	apf::destroyMeshElement(me);
	me = NULL;
	ke2 = stiff2.ke;

	std::size_t nrows1, ncols1, nrows2, ncols2;
	nrows1 = ke1.getRows();
	ncols1 = ke1.getColumns();
	nrows2 = ke2.getRows();
	ncols2 = ke2.getColumns();

	assert(nrows1 == nrows2);
	assert(ncols1 == ncols2);

	ke3.setSize(nrows1, ncols1);

	for(std::size_t ii = 0; ii < nrows1; ++ii) {
		for(std::size_t jj = 0; jj < ncols1; ++jj) {
			float ratio = ke1(ii,jj) / ke2(ii, jj);
			if(fabs(ratio - 1.0) > 1e-6) {
				ke3(ii, jj) = 1;
			} else {
				ke3(ii, jj) = 8;
			}
			EXPECT_NEAR(ke1(ii, jj), ke2(ii, jj), 5e-15);
		}
	}
}

TEST_F(LocalStiffnessMatrixTest, GetMatrix) {
	apf::Mesh2* two_mesh = NULL;
	this->mesh_builder->build2DRectQuadMesh(two_mesh, 1, 2, 0.0, 0.0, 2.0, 1.0);
	apf::changeMeshShape(two_mesh, apf::getLagrange(2));
	
	this->field1 = createField(two_mesh, "dummy", apf::VECTOR, two_mesh->getShape());
	apf::zeroField(this->field1); /*must zero the field to force writing of data*/
	EXPECT_TRUE(two_mesh != NULL);


	if(two_mesh != NULL) {
		two_mesh->destroyNative();
		apf::destroyMesh(two_mesh);
	}
}

/*=======================================================================*/
class UserFieldTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
		mesh_builder = new MeshBuilder();
		mesh = NULL;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		delete mesh_builder;	
	}

};

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
			apf::destroyMeshElement(me);
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

TEST_F(UserFieldTest, ScratchPad) {
	/*Creats a user defined field aboive in class MyUserFunction 
	* and checks the gradient of that field at several randomly
	* chosen points within the boundary of that element*/
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
	/*cannot destroy mesh element until field element is finished*/
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

	/*clean up the field element*/
	apf::destroyElement(f_elm);
	apf::destroyMeshElement(me);
	me = NULL;
	apf::destroyField(test_f);

	apf::writeASCIIVtkFiles("batman_elm", m);

	m->destroyNative();
	apf::destroyMesh(m);
}

TEST_F(StiffnessTest, QuadraticShapeFunctionCheck) {



}
