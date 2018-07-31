#ifndef AddOns_Recola_Recola_Born_H
#define AddOns_Recola_Recola_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace Recola {

  class Recola_Born: public PHASIC::Tree_ME2_Base {

    int m_recola_id, m_amptype;
    double m_symfac;

  public:

    Recola_Born(const PHASIC::External_ME_Args& args,
		unsigned int recola_id, int amptype);
    
    ~Recola_Born() {}
    
    double Calc(const ATOOLS::Vec4D_Vector& momenta);
    
  };
  
}

#endif
