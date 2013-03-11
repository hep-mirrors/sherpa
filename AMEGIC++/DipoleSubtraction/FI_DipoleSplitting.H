#ifndef AMEGIC_DipoleSubtraction_FI_DipoleSplitting_H
#define AMEGIC_DipoleSubtraction_FI_DipoleSplitting_H

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

  class FI_DipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_zi,m_zj;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    FI_DipoleSplitting(int ft, int m,int i,int j,int k) 
    {m_ft=ft; m_type=dpt::f_i; m_m=m; m_i=i;m_j=j;m_k=k;}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return 1.-m_xijk;}
    void CalcDiPolarizations();
  };

  class FI_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_zi,m_zj;
    double m_mi,m_mj,m_mij;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    FI_MassiveDipoleSplitting(int ft, int m,int i,int j,int k, double mi, double mj, double mij) 
    {m_ft=ft; m_type=dpt::f_im; m_m=m; m_i=i;m_j=j;m_k=k;
    m_mi=mi*mi; m_mj=mj*mj; m_mij=mij*mij; }
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return 1.-m_xijk;}
    void CalcDiPolarizations();
  };
}
#endif
