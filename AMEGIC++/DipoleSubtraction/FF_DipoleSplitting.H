#ifndef AMEGIC_DipoleSubtraction_FF_DipoleSplitting_H
#define AMEGIC_DipoleSubtraction_FF_DipoleSplitting_H

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

  class FF_DipoleSplitting : public DipoleSplitting_Base {
    double m_yijk;
    double m_zi,m_zj;
    double m_sff;
//     ATOOLS::Vec4D m_pt;
  public:
    FF_DipoleSplitting(int ft, int m,int i,int j,int k) 
    {m_ft=ft; m_type=dpt::f_f; m_m=m; m_i=i;m_j=j;m_k=k;}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetF();
    double GetDPSF()   {return m_yijk;}
    void CalcDiPolarizations();
  };

  class FF_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_yijk;
    double m_zi,m_zj;
    double m_zim,m_zjm;
    double m_sff;
    double m_mi,m_mj,m_mk,m_mij,m_Q2;
    double m_vijk,m_zpm;
    ATOOLS::Vec4D m_Q;
  public:
    FF_MassiveDipoleSplitting(int ft, int m,int i,int j,int k, double mi, double mj, double mk, double mij) 
    { m_ft=ft; m_type=dpt::f_fm; m_m=m; m_i=i;m_j=j;m_k=k;
      m_mi=mi*mi; m_mj=mj*mj; m_mk=mk*mk; m_mij=mij*mij; }
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetF();
    double GetDPSF()   {return m_yijk;}
    void CalcDiPolarizations();
  };
}
#endif
