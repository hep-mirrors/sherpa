#ifndef AMISIC_Perturbative_QCD_Processes_H
#define AMISIC_Perturbative_QCD_Processes_H

#include "AMISIC++/Perturbative/MI_Process.H"
#include "ATOOLS/Phys/Flavour.H"

namespace AMISIC {
  class qg_qgamma: public XS_Base {
  public:
    qg_qgamma();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };
  
  class qqbar_ggamma: public XS_Base {
  public:
    qqbar_ggamma();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };
}

#endif
