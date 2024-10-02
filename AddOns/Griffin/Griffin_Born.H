#ifndef AddOns_Griffin_Griffin_Born_H
#define AddOns_Griffin_Griffin_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/Griffin/Griffin_Interface.H"

namespace MODEL {
    class Model_Base;
}

namespace Griffin {

    class Griffin_Born: public PHASIC::Tree_ME2_Base {

    int m_Griffin_id;
    int m_amptype;
    double m_symfac;

    public:

    Griffin_Born(const PHASIC::External_ME_Args& args);
                              
    ~Griffin_Born() {}
                                  
    double Calc(const ATOOLS::Vec4D_Vector& momenta);
                                      
    };
      
}

#endif

