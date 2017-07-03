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
  public:
    IF_DipoleSplitting(ATOOLS::sbt::subtype st,
                       ATOOLS::spt::splittingtype ft,
                       int m,int i,int j,int k) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::i_f,m,i,j,k) {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_uj;}
    void CalcDiPolarizations();
  };

  class IF_MassiveDipoleSplitting : public DipoleSplitting_Base {
    double m_xijk;
    double m_uj,m_uk;
  public:
    IF_MassiveDipoleSplitting(ATOOLS::sbt::subtype st,
                              ATOOLS::spt::splittingtype ft,
                              int m,int i,int j,int k)  :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::i_fm,m,i,j,k) {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_uj;}
    void CalcDiPolarizations();
  };
}
#endif
