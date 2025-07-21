#ifndef HADRON_RESCATTERING_XSecs_MesonMeson_H
#define HADRON_RESCATTERING_XSecs_MesonMeson_H

#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include <map>
#include <list>

namespace HADRON_RESCATTERING {
  class MesonMeson {
    std::map<HR_Res, HR_resonances>  m_resonances;    

    void InitialiseResonances();
    
    bool m_test;
    void Tests();
  public:
    MesonMeson();
    ~MesonMeson();

    double XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
		 const double & s);
    double XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
		const double & s);
  };
}

#endif

