#ifndef CSSHOWER_Sudakov_CS_SingleSudakov_QCD_H
#define CSSHOWER_Sudakov_CS_SingleSudakov_QCD_H

#include "CSSHOWER++/Sudakov/CS_SingleSudakov_Base.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

namespace CSSHOWER {
  //the FF Sudakovs 
  class DeltaFF_q_qg : public CS_SingleSudakov_Base {
    double m_mi2;
    ATOOLS::Gauss_Integrator m_gauss;
  
  public:
    DeltaFF_q_qg(const ATOOLS::Flavour &);
    double operator()(double kt2prod,double kt2dec);
    double operator()(double kt2dec);
    double Delta(double kt2prod,double kt2dec);
  };

  class DeltaFF_g_qq : public CS_SingleSudakov_Base {
    double m_mj2;
    ATOOLS::Gauss_Integrator m_gauss;
  public:
    DeltaFF_g_qq(const ATOOLS::Flavour &);
    double operator()(double kt2prod,double kt2dec);
    double operator()(double kt2dec);
    double Delta(double kt2prod,double kt2dec);
  };

  class DeltaFF_g_gg : public CS_SingleSudakov_Base {
    ATOOLS::Gauss_Integrator m_gauss;
  public:
    DeltaFF_g_gg();
    double operator()(double kt2prod,double kt2dec);
    double operator()(double kt2dec);
    double Delta(double kt2prod,double kt2dec);
  };

}

#endif
