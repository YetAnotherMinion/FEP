#include <iostream>

#include <apfMesh.h>

#include "ForceContributor2D.h"

ForceContributor2D::ForceContributor2D(
    apf::Field* f,
    uint32_t integrate_order,
    apf::Vector3(*fnc)(apf::Vector3 const& p)) : apf::Integrator(integrate_order), field(f)
{
    this->ndims = apf::getMesh(f)->getDimension();
    this->fnc = fnc;
}

void ForceContributor2D::inElement(apf::MeshElement* me)
{
    /*create the Field Element from the MeshElement*/
    this->field_element = apf::createElement(this->field, me);
    /*determine the size of the force matrices*/
    this->nnodes = apf::countNodes(this->field_element);
    this->ndofs = this->ndims * apf::countNodes(this->field_element);
    this->fe.reserve(this->ndofs);
    /*zero the internal matricies*/
    this->fe.clear();
    fe.resize(this->ndofs, 0.0);

}

void ForceContributor2D::outElement()
{
    std::cout << "==========================" << std::endl;
    for(std::size_t ii = 0; ii < this->fe.size(); ++ii) {
        std::cout << ii << " : " << this->fe[ii] << std::endl;
    }
    /*perform clean up and destroy specific field element*/
    apf::destroyElement(this->field_element);
}

void ForceContributor2D::atPoint(apf::Vector3 const& p, double w, double dV)
{
    apf::NewArray<double> shape_val;
    apf::getShapeValues(this->field_element, p, shape_val);

    apf::Vector3 x;
    apf::MeshElement* me = apf::getMeshElement(this->field_element);
    apf::mapLocalToGlobal(me, p, x);
    /*compute the contributions of each shape function*/
    apf::Vector3 tmp = this->fnc(x);
    for(uint32_t ii = 0; ii < this->nnodes; ++ii){
        for(uint32_t jj = 0; jj < this->ndims; ++jj){
            this->fe[(ii*this->ndims)+jj] += (shape_val[ii] * tmp[jj] * w * dV);
        }
    }
}
