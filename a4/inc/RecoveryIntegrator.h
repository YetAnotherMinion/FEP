#ifndef FEP_RECOVERY_INTEGRATOR_H
#define FEP_RECOVERY_INTEGRATOR_H 1

#include <stdint.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>

/*template the size of the D to enable this function to*/
template< uint32_t N >
class RecoveryIntegrator : public apf::Integrator
{
public:
    /* Takes a function that computes the D for a point*/
	RecoveryIntegrator(
        apf::Field *f,
        apf::Numbering *,
        const apf::Matrix< N, N >& (*)(apf::Vector3),
        uint32_t integrate_order,
        uint32_t n_global_dofs);

	void inElement(apf::MeshElement *me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dV);
    apf::DynamicMatrix ke;
private:
    apf::Field* field;
    apf::Element* field_element;
    uint32_t ndofs; /*gets updated per element processed*/
    uint32_t nnodes;
    uint32_t ndims;
    uint32_t n_global_dofs;
    const apf::Matrix< N, N >& (*ansiotropy)(apf::Vector3);
    std::vector<double> stress;
    std::vector<double> strain;
    std::vector<double> energy; /*will hold the strain energy per element*/
    
};

#endif
