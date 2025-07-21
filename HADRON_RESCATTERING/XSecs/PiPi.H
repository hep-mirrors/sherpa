#ifndef HADRON_RESCATTERING_XSecs_PiPi_H
#define HADRON_RESCATTERING_XSecs_PiPi_H

#include "ATOOLS/Phys/Flavour.H"

//////////////////////////////////////////////////////////////////////////////
//
// All total cross sections in mb.
// For details regarding parametrizations etc., cf the C-file.
//
//////////////////////////////////////////////////////////////////////////////

namespace HADRON_RESCATTERING {  
  class PiPi {
    double m_offset;
  public:
    PiPi();
    ~PiPi() {}
    double XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
		 const double & s);
    double XSel(const double & s);
  };
}

#endif
