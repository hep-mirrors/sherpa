#ifndef AMEGIC_DipoleSubtraction_IF_DipoleSplitting_H
#define AMEGIC_DipoleSubtraction_IF_DipoleSplitting_H

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

  class IF_DipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_uj,m_uk;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    IF_DipoleSplitting(int ft, int m,int i,int j,int k) 
    {m_ft=ft; m_type=dpt::i_f; m_m=m; m_i=i;m_j=j;m_k=k;}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetF();
    double GetDPSF()   {return m_uj;}
    void CalcDiPolarizations();
  };

  class IF_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_uj,m_uk;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    IF_MassiveDipoleSplitting(int ft, int m,int i,int j,int k) 
    {m_ft=ft; m_type=dpt::i_fm; m_m=m; m_i=i;m_j=j;m_k=k;}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetF();
    double GetDPSF()   {return m_uj;}
    void CalcDiPolarizations();
  };
}
#endif
