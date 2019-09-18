#ifndef AddOns_Recola_Recola_Born_H
#define AddOns_Recola_Recola_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/Recola/Recola_Interface.H"

namespace MODEL {
  class Model_Base;
}

namespace Recola {

  class Recola_Born: public PHASIC::Tree_ME2_Base {

    int m_recola_id;
    int m_ewscheme;
    double m_symfac;
    size_t m_oqcd;

  public:

    Recola_Born(const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs,
		unsigned int recola_id);
    
    ~Recola_Born() {}
    
    double Calc(const ATOOLS::Vec4D_Vector& momenta);
    
  };
  
}

#endif
