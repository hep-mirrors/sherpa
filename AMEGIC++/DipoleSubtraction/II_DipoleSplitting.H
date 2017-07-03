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
  public:
    II_DipoleSplitting(ATOOLS::sbt::subtype st,
                       ATOOLS::spt::splittingtype ft,
                       int m,int i,int j,int k) :
      DipoleSplitting_Base(st,ft,ATOOLS::dpt::i_i,m,i,j,k) {}
    
    void SetMomenta(const ATOOLS::Vec4D*);
    double GetValue();
    double GetDPSF()   {return m_vi;}
    void CalcDiPolarizations();
  };
}
#endif
