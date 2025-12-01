#ifndef HADRON_RESCATTERING_XSecs_BaryonMeson_H
#define HADRON_RESCATTERING_XSecs_BaryonMeson_H


#include "HADRON_RESCATTERING/XSecs/NNfits.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "ATOOLS/Phys/Flavour.H"
#include "HADRON_RESCATTERING/XSecs/NPion.H"


namespace HADRON_RESCATTERING {
  class BaryonMeson {
    NNfits  m_BM;
    NPion m_NP;
    
    HPR1R2  m_hpr1r2;
    double  m_s, m_mA, m_mA2, m_mB, m_mB2;
    
    double xstot(long int & A, long int & B, const double & plab);
    double xsel(long int & A, long int & B, const double & plab);
    
    bool m_test = false;
    void Tests();
  public:
    BaryonMeson();
    ~BaryonMeson();

    double XStot(const ATOOLS::Flavour & A, const ATOOLS::Flavour & B, const double & s);
    double XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B, const double & s);
  };
}








#endif