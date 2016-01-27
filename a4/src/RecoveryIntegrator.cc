#include <apfMesh.h>

#include "RecoveryIntegrator.h"

template<unsigned int N> RecoveryIntegrator<N>::RecoveryIntegrator(
    apf::Field* f,
    apf::Numbering* nodeNums,
    const apf::Matrix< N,N > & (*stiffness_function)(apf::Vector3),
    uint32_t integrate_order,
    uint32_t n_global_dofs) : apf::Integrator(integrate_order), field(f),
                     ansiotropy(stiffness_function), n_global_dofs(n_global_dofs)
{
    stress.resize(n_global_dofs);
    strain.resize(n_global_dofs);
    ndims = apf::getMesh(f)->getDimension();
}

template<unsigned int N> void RecoveryIntegrator<N>::inElement(apf::MeshElement* me)
{
    /*create the Field Element from the MeshElement*/
    this->field_element = apf::createElement(this->field, me);
    /*determine the size of the force matrices*/
    this->nnodes = apf::countNodes(this->field_element);
    this->ndofs = this->ndims * apf::countNodes(this->field_element);
}

template<unsigned int N> void RecoveryIntegrator<N>::outElement()
{
    /*perform clean up and destroy specific field element*/
    apf::destroyElement(this->field_element);
}


template<unsigned int N> void RecoveryIntegrator<N>::atPoint(apf::Vector3 const& p, double w, double dV)
{
    // below two lines do not appear to be used [Jan 14 2016]
    // apf::NewArray<double> shape_val;
    // apf::getShapeValues(this->field_element, p, shape_val);

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

        }
    }

}
