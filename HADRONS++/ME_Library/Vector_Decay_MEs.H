#ifndef HADRONS_ME_Library_Vector_Decay_MEs_H
#define HADRONS_ME_Library_Vector_Decay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {
  class V_FF : public HD_ME_Base {
  public:
    V_FF(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class V_PP : public HD_ME_Base {
  public:
    V_PP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class V_VP : public HD_ME_Base {
  public:
    V_VP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class V_AP : public HD_ME_Base {
  public:
    V_AP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class V_VS : public HD_ME_Base {
  public:
    V_VS(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class V_VV : public HD_ME_Base {
  public:
    V_VV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class A_VP : public HD_ME_Base {
  public:
    A_VP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };


  class V_PFF : public HD_ME_Base {
  public:
    V_PFF(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };
  
  class V_PPP : public HD_ME_Base {
  public:
    V_PPP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };
  
  /*!
    \class V_FF
    \brief For decays \f$V\to\ell\bar\ell\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_\mu\bar u_\ell\gamma^\mu u_{\bar\ell}
    \f]
  */
  /*!
    \class V_PP
    \brief For decays \f$V\to\pi\pi\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_\mu(p_1-p_2)^\mu
    \f]
  */
  /*!
    \class V_VP
    \brief For decays \f$V\to V'\pi\f$ and similar
    
    \f[
    {\cal M} = g\epsilon^{\mu\nu\rho\sigma}
    \epsilon_\mu p_\nu\epsilon'_\rho p'_\sigma
    \f]
  */
  /*!
    \class V_AP
    \brief For decays \f$V\to A\pi\f$ and similar
    
    \f[
    {\cal M} = g\left[(p_V\cdot p_A)(\epsilon_V\cdot\epsilon_A)-
                      (\epsilon_V\cdot p_A)(p_V\cdot\epsilon_A)\right]
    \f]
  */
  /*!
    \class V_VS
    \brief For decays \f$V\to V'S\f$ and similar
    
    \f[
    {\cal M} = g\epsilon^\mu\left[(p\cdot p')\epsilon'_\mu-
                                  (p\cdot \epsilon')p'_\mu\right]

    \f]
  */
  /*!
    \class V_VV
    \brief For decays \f$V_0\to V_1V_2\f$ and similar
    
    \f[
    {\cal M} = \frac{g}{M_0^2}\left[
    \left(p_{0,\mu}\epsilon_{0,\nu}-p_{0,\nu}\epsilon_{0,\mu}\right)
    \left(p_1^\mu\epsilon_1^\rho-p_1^\rho\epsilon_1^\mu\right)
    \left(p_2^\nu\epsilon_{2,\rho}-p_{2,\rho}\epsilon_2^\nu\right) +
    \left(\mu\leftrightarrow\nu\right)
    \right]
    \f]
  */
  /*!
    \class A_VP
    \brief For decays \f$a_1\to \rho\pi\f$ and similar
    
    \f[
    {\cal M} = g\epsilon_{\mu,A}
    \left[(p_V\cdot p_A)\epsilon_V^\mu-(\epsilon_V\cdot p_A)p_V^\mu\right]
    \f]
  */
  /*!
    \class V_PFF
    \brief For decays \f$\rho\to\pi\ell\bar\ell\f$ and similar
    
    \f[
    {\cal M} = \frac{g}{q^2_{\ell\bar\ell}}\epsilon_{\mu\nu\rho\sigma}
    \epsilon_V^\mu p_V^\nu q_{\ell\bar\ell}^\rho
    \bar u_\ell\gamma^\sigma u_{\bar\ell}
    \f]
  */
  /*!
    \class V_PPP
    \brief For decays \f$\omega\to\pi\pi\pi\f$ and similar
    
    \f[
    \f]
    cf.\ hep-ex/9902022, hep-ex/9908024 and hep-ex/0303016.
  */
};

#endif
