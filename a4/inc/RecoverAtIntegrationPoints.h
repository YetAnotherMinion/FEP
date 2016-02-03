#ifndef FEP_RECOVERY_INTEGRATOR_H
#define FEP_RECOVERY_INTEGRATOR_H 1

#include <stdint.h>
#include <utility>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>
#include <apfMesh.h>

/*template the size of the D to enable this function to work for
* 2D and 3D problems*/
template< uint32_t N, uint32_t M_COMPONENTS >
class RecoverAtIntegrationPoints : public apf::Integrator
{
public:
	RecoverAtIntegrationPoints(
        apf::Field *f,
        const apf::Matrix< N, N >,
        std::vector< apf::Vector<M_COMPONENTS> > &,
        uint32_t integrate_order);

	void inElement(apf::MeshElement *me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dV);
    void generate_B(apf::Element*, apf::NewArray < apf::Matrix< N, M_COMPONENTS > > & B, apf::Vector3 const & point);

    /*the first member is the location in space, the second member*/
    std::vector< std::pair<apf::Vector3, apf::Vector<N> > > strain;
    std::vector< std::pair<apf::Vector3, apf::Vector<N> > > stress;
    std::vector<double> energy; /*will hold the strain energy per element*/
private:
    apf::Field* field;
    apf::Element* field_element;
    apf::Matrix< N, N> D;
    uint32_t ndofs; /*gets updated per element processed*/
    uint32_t nnodes;
    uint32_t ndims;
    std::vector< apf::Vector<M_COMPONENTS> > & local_displacements;
};
/*----------Implementation----------------*/

template<uint32_t N, uint32_t M> RecoverAtIntegrationPoints<N, M>::RecoverAtIntegrationPoints(
    apf::Field* f,
    const apf::Matrix< N,N > D,
    std::vector<apf::Vector<M> > & local_disp,
    uint32_t integrate_order) : apf::Integrator(integrate_order), field(f), D(D),
                     local_displacements(local_disp)
{
    /*we could reserve space in the strain and stress vectors here
    * if we wanted to find the numer*/
    ndims = apf::getMesh(f)->getDimension();
}

template<uint32_t N, uint32_t M> void RecoverAtIntegrationPoints<N, M>::inElement(apf::MeshElement* me)
{
    /*create the Field Element from the MeshElement*/
    this->field_element = apf::createElement(this->field, me);
    /*determine the size of the force matrices*/
    this->nnodes = apf::countNodes(this->field_element);
    this->ndofs = this->ndims * apf::countNodes(this->field_element);
}

template<uint32_t N, uint32_t M> void RecoverAtIntegrationPoints<N, M>::outElement()
{
    /*perform clean up and destroy specific field element*/
    apf::destroyElement(this->field_element);
}


template<uint32_t N, uint32_t M> void RecoverAtIntegrationPoints<N, M>::atPoint(apf::Vector3 const& p, double w, double dV)
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

    apf::NewArray< apf::Matrix< N, M > > B(this->nnodes);
    this->generate_B(field_element, B, p);

    /*compute the strain at this point*/
    apf::Vector<N> strain_at_point;
    for(uint32_t ii = 0; ii < this->nnodes; ++ii) {
        apf::Vector<N> foo;

    }
    this->strain.push_back(std::make_pair(p, strain_at_point));
    apf::Vector<N> stress_at_point;
    stress_at_point = D * strain_at_point;
    this->stress.push_back(std::make_pair(p, stress_at_point));
}

/*throw exception if we have not specialized a template for a specific
* dimensionality problem*/
template<uint32_t N, uint32_t M> void RecoverAtIntegrationPoints<N, M>::generate_B(apf::Element* field_element, apf::NewArray < apf::Matrix<N,M> > & B, apf::Vector3 const & p)
{
    throw std::runtime_error("Dimensionality not implemented");
}

// /*template specifications for 2D linear elastic problem*/
// template<3, 2> void RecoverAtIntegrationPoints<N, M>::generate_B(apf::Element* field_element, apf::NewArray < apf::Matrix<N, M> > & B, apf::Vector3 const & point)
// {
//     apf::NewArray<apf::Vector3> gradShape;
//     apf::getShapeGrads(field_element, point, gradShape);
//     /*construct each of the nnodes shape function matricies*/
//     uint32_t ii, jj;
//     for(ii = 0; ii < this->nnodes; ++ii) {
//         B[ii][0][0] = gradShape[ii][0];
//         B[ii][0][1] = 0.0;
//         B[ii][1][0] = 0.0;
//         B[ii][1][1] = gradShape[ii][1];
//         B[ii][2][0] = gradShape[ii][1];
//         B[ii][2][1] = gradShape[ii][0];
//     }
// }

// /*template specification for 3D linear elastic problem*/
// template<6, 3> void RecoverAtIntegrationPoints<N, M>::generate_B(apf::Element* field_element, apf::NewArray < apf::Matrix<N, M> > & B, apf::Vector3 const & p)
// {
    
// }

#endif
