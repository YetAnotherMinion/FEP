#include <apfMesh.h>

#include "StiffnessContributor2D.h"

StiffnessContributor2D::StiffnessContributor2D(
	apf::Field* f,
	apf::Matrix< 3,3 > D_input,
	uint32_t integrate_order) : apf::Integrator(integrate_order), field(f), D(D_input)
{
  ndims = apf::getMesh(f)->getDimension();
}

void StiffnessContributor2D::inElement(apf::MeshElement* me)
{
	/*create the Field Element from the MeshElement*/
	this->field_element = apf::createElement(this->field, me);
	/*determine the size of the force matrices*/
	this->nnodes = apf::countNodes(this->field_element);
	this->ndofs = this->ndims * apf::countNodes(this->field_element);
	this->ke.setSize(this->ndofs, this->ndofs);
	/*zero the internal matricies*/
	this->ke.zero();
}

void StiffnessContributor2D::outElement()
{
	// std::cout << this->ke << std::endl;
	for(uint32_t ii = 0; ii < this->ndofs; ++ii) {
		for(uint32_t jj = ii; jj < this->ndofs; ++jj) {
			if(ii != jj) {
				this->ke(jj,ii) = this->ke(ii, jj);
			}
		}
	}
	/*perform clean up and destroy specific field element*/
	apf::destroyElement(this->field_element);
}

void StiffnessContributor2D::atPoint(apf::Vector3 const& p, double w, double dV)
{
	apf::NewArray<double> shape_val;
	apf::getShapeValues(this->field_element, p, shape_val);

	apf::NewArray<apf::Vector3> gradShape;
	apf::getShapeGrads(this->field_element, p, gradShape);

	/*the below three lines do not appear to be used for anything [Dec 12 2015]*/
	// apf::Vector3 x;
	// apf::MeshElement* me = apf::getMeshElement(this->field_element);
	// apf::mapLocalToGlobal(me, p, x);

	apf::NewArray< apf::Matrix< 3,2 > > B(this->nnodes);
	/*construct each of the nnodes shape function matricies*/
	uint32_t ii, jj;
	for(ii = 0; ii < this->nnodes; ++ii) {
		B[ii][0][0] = gradShape[ii][0];
		B[ii][0][1] = 0.0;
		B[ii][1][0] = 0.0;
		B[ii][1][1] = gradShape[ii][1];
		B[ii][2][0] = gradShape[ii][1];
		B[ii][2][1] = gradShape[ii][0];
	}

	/*assemble all blocks since there are inconsistent results when trying to assemble
	* only above main diagonal*/
	for(ii = 0; ii < this->nnodes; ++ii) {
		for(jj = ii; jj < this->nnodes; ++jj) {
			apf::Matrix< 2,2 > nodal_submatrix = (transpose(B[ii]) * (this->D) * B[jj]);
			nodal_submatrix = nodal_submatrix * w * dV;
			/*add the contribution to the element stiffness matrix,
			* we unroll the element access loop since it is so small*/
			this->ke(2*ii,(2*jj)) += nodal_submatrix[0][0];
			this->ke(2*ii,(2*jj) + 1) += nodal_submatrix[0][1];
			if(ii != jj) {
				this->ke((2*ii) + 1,2*jj) += nodal_submatrix[1][0];
			}
			this->ke((2*ii) + 1,2*jj + 1) += nodal_submatrix[1][1];
		}
	}

}
