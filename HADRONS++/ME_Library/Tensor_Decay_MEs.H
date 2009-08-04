#ifndef HADRONS_ME_Library_Tensor_Decay_MEs_H
#define HADRONS_ME_Library_Tensor_Decay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {
  class T_PP : public HD_ME_Base {
  public:
    T_PP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class T_VP : public HD_ME_Base {
  public:
    T_VP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class T_VV : public HD_ME_Base {
  public:
    T_VV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  /*!
    \class T_PP
    \brief For decays \f$T\to\pi\pi\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_{\mu\nu}p_1\mu p_2^\nu
    \f]
  */
  /*!
    \class T_VP
    \brief For decays \f$T\to V\pi\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_T^{\mu\nu}p_{P,\mu}\epsilon_{\nu\rho\sigma\kappa}
    p_V^\nu\epsilon_V^\rho p_P^\kappa
    \f]
  */
  /*!
    \class T_VV
    \brief For decays \f$T\to VV'\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_T^{\mu\nu}
    \left(p_\mu\epsilon_\rho-p_\rho\epsilon_\mu\right)
    \left(p'_\nu{\epsilon'}^\rho-{p'}^\rho\epsilon'_\nu\right)
    -\left(\mu\leftrightarrow\nu\right)
    \f]
  */
};

#endif
