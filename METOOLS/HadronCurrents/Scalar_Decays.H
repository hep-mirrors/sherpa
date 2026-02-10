#ifndef METOOLS_HadronCurrents_Scalar_Decays_H
#define METOOLS_HadronCurrents_Scalar_Decays_H

#include "METOOLS/HadronCurrents/Resonance_Base.H"

namespace METOOLS {
  ////////////////////////////////////////////////////////////////////
  //
  // (Pseudo)Scalar -> (Pseudo)Scalar + (Pseudo)Scalar
  //
  ////////////////////////////////////////////////////////////////////
  class S_PP: public Partial_Width_Base {
    double m_lambda2;
    const double ME2(const double & s) const;
    const double PS_2(const double & s) const;

    const double Calculate(const double & s);
  public:
    S_PP(const ATOOLS::Flavour & inflav,
	 const std::vector<ATOOLS::Flavour> & outflavs,
	 const double & BR);    
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange) {};
  };
}

#endif
