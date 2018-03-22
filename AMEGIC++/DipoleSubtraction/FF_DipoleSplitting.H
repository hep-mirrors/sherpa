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
    double m_vi_tilde, m_zain, m_viab, m_xain; // for ID-dipoles
    ATOOLS::Vec4D m_n;                         // for ID-dipoles
  public:
    FF_DipoleSplitting(ATOOLS::sbt::subtype st,
                       ATOOLS::spt::splittingtype ft,
                       int m,int i,int j,int k) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::f_f,m,i,j,k),
      m_yijk(0.), m_zi(0.), m_zj(0.) 
    {}

    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_yijk;}
    void CalcDiPolarizations();
  };

  class FF_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_yijk, m_yp;
    double m_zi,m_zj;
    double m_zim,m_zjm;
    double m_mi2,m_mj2,m_mk2,m_mij2,m_Q2;
    double m_vijk,m_zpm;
    ATOOLS::Vec4D m_Q;
  public:
    FF_MassiveDipoleSplitting(ATOOLS::sbt::subtype st,
                              ATOOLS::spt::splittingtype ft,
                              int m,int i,int j,int k,
                              double mi, double mj, double mk, double mij) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::f_fm,m,i,j,k),
      m_yijk(0.), m_yp(0.), m_zi(0.), m_zj(0.), m_zim(0.), m_zjm(0.),
      m_mi2(mi*mi), m_mj2(mj*mj), m_mk2(mk*mk), m_mij2(mij*mij),
      m_Q2(0.), m_vijk(0.), m_zpm(0.), m_Q(ATOOLS::Vec4D(0.,0.,0.,0.))
    {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_yijk;}
    void CalcDiPolarizations();
  };
}
#endif
