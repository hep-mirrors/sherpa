#ifndef METOOLS_HadronCurrents_Scalar_Decays_H
#define METOOLS_HadronCurrents_Scalar_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

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

  ////////////////////////////////////////////////////////////////////
  //
  // Generic placeholder for a multi-body channel with no dedicated
  // phase-space/Dalitz treatment implemented yet: contributes a
  // CONSTANT partial width above the sum-of-masses threshold (using
  // the base class's own threshold theta-function, which works for
  // any number of outflavs, not just two), and exactly zero below it.
  // This is NOT a real dynamical model of the channel's s-dependence -
  // it exists purely to keep the total width budget honest (e.g. for
  // a0(1450)'s dominant omega pi pi channel, ~67% of the width, for
  // which building a proper 3-body S->VPP Dalitz integration is out of
  // scope for now). Flag prominently wherever used.
  //
  ////////////////////////////////////////////////////////////////////
  class Constant_PW: public Partial_Width_Base {
    const double Calculate(const double & s) { return m_prefactor; }
  public:
    Constant_PW(const ATOOLS::Flavour & inflav,
		const std::vector<ATOOLS::Flavour> & outflavs,
		const double & BR) :
      Partial_Width_Base(inflav,outflavs,BR)
    {
      FixPrefactor();
    }
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange) {};
  };
}

#endif
