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


  class II_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_xijk, m_vi, m_sff, m_xab, m_vab;
    double m_ma, m_mb, m_mk;
    double m_xmin;
    int i_i,i_j,i_k;
    int m_test;
  public:
    II_MassiveDipoleSplitting(int ft,
			      int m,int i,int j,int k,
			      double mi, double mj, double mk) 
    {
      m_ft=ft; m_type=dpt::i_im; m_m=m; m_i=i; m_j=j; m_k=k;
      i_i=i;i_j=j;i_k=k; m_ma=mi*mi;m_mb=mk*mk; m_mk=mj*mj;
    }
    void   SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()        {return m_xab;}
    double GetXMin()        {return m_xmin;}
    void   CalcDiPolarizations();
    inline double lambda(double a, double b, double c) {return ATOOLS::sqr(a-(b+c))-4.*b*c ;}
  };
}
#endif
