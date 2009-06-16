#ifndef HADRONS_ME_Library_Phi_Decay_MEs_H
#define HADRONS_ME_Library_Phi_Decay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {
  class Phi_PPGamma : public HD_ME_Base {
    public:
      Phi_PPGamma(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
	HD_ME_Base(flavs,n,indices,name) {};
      void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
      void SetModelParameters(GeneralModel);
  };

  class Phi_VP : public HD_ME_Base {
    public:
      Phi_VP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
	HD_ME_Base(flavs,n,indices,name) {};
      void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
      void SetModelParameters(GeneralModel);
  };

  /*!
    \class Phi_PPGamma
    \brief For decays \f$\phi\to\pi\pi\gamma\f$, \f$\phi\to\eta\pi\gamma\f$

    \f[
    \f]
    cf.\ hep-ph/0404076, hep-ph/0306249, hep-ph/0212153.
   */
  /*!
    \class Phi_VP
    \brief For decays \f$\phi\to\omega\pi\f$

    \f[
    \f]
    cf.\ hep-ph/0404072.
   */

  // btw: D_s->phi form factor in hep-ph/0308259

};

#endif
