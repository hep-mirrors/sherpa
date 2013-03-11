#ifndef AMEGIC_DipoleSubtraction_II_DipoleSplitting_H
#define AMEGIC_DipoleSubtraction_II_DipoleSplitting_H

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

  class II_DipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_vi;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    II_DipoleSplitting(int ft, int m,int i,int j,int k) 
    {m_ft=ft; m_type=dpt::i_i; m_m=m; m_i=i;m_j=j;m_k=k;}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_vi;}
    void CalcDiPolarizations();
  };
}
#endif
