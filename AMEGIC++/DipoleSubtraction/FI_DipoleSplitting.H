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
  public:
    FI_DipoleSplitting(ATOOLS::sbt::subtype st,
                       ATOOLS::spt::splittingtype ft,
                       int m,int i,int j,int k) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::f_i,m,i,j,k) {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return 1.-m_xijk;}
    void CalcDiPolarizations();
  };

  class FI_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_zi,m_zj;
    double m_mi,m_mj,m_mij;
  public:
    FI_MassiveDipoleSplitting(ATOOLS::sbt::subtype st,
                              ATOOLS::spt::splittingtype ft,
                              int m,int i,int j,int k,
                              double mi, double mj, double mij) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::f_im,m,i,j,k),
      m_mi(mi*mi), m_mj(mj*mj), m_mij(mij*mij) {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return 1.-m_xijk;}
    void CalcDiPolarizations();
  };
}
#endif
