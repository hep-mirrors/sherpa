#ifndef HADRONS_ME_Library_Baryon_SimpleDecay_MEs_H
#define HADRONS_ME_Library_Baryon_SimpleDecay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {
  class D_Radiative_E1 : public HD_ME_Base {
  public:
    D_Radiative_E1(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class D_Radiative_M1 : public HD_ME_Base {
  public:
    D_Radiative_M1(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class R_Radiative_E1 : public HD_ME_Base {
  public:
    R_Radiative_E1(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class R_Radiative_M1 : public HD_ME_Base {
  public:
    R_Radiative_M1(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };


  class D_DP : public HD_ME_Base {
  public:
    D_DP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class D_DV : public HD_ME_Base {
  public:
    D_DV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class D_RP : public HD_ME_Base {
  public:
    D_RP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class D_RV : public HD_ME_Base {
  public:
    D_RV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class R_DP : public HD_ME_Base {
  public:
    R_DP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class R_DV : public HD_ME_Base {
  public:
    R_DV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  class R_RP : public HD_ME_Base {
  public:
    R_RP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,HELICITIES::Spin_Amplitudes *,int);
    void SetModelParameters(GeneralModel);
  };

  /*!
    \class D_Radiative_E1
    \brief For radiative decays \f$\frac12\to \frac{1}{2'}\gamma\f$ through
    an electric transition.

    \f[
    {\cal M} = A_{E1}\bar u_{1/2'}\gamma^\mu\gamma_5u_{1/2}F^{\mu\nu}p_{0\mu}
    \f]
  */

  /*!
    \class D_Radiative_M1
    \brief For radiative decays \f$\frac12\to \frac{1}{2'}\gamma\f$ through
    a magnetic transition.

    \f[
    {\cal M} = A_{M1}\bar u_{1/2'}\sigma^{\mu\nu}\bar u_{1/2}F^{\mu\nu}\f]
    \f]

    Here, the Gordon-identity can/must be used to rephrase this in terms
    of known expressions for the spinor line.
  */

  /*!
    \class R_Radiative_E1
    \brief For radiative decays \f$\frac32\to \frac12\gamma\f$ through
    an electric transition.

    \f[
    {\cal M} = A_{E1}\bar u_{1/2}\bar u_{\nu, 3/2} F^{\mu\nu}p_{0\mu}
    \f]
  */

  /*!
    \class R_Radiative_M1
    \brief For radiative decays \f$\frac32\to\frac12\gamma\f$ through
    a magnetic transition.

    \f[
    {\cal M} = A_{M1}\bar u_{1/2}\gamma_\mu\gamma_5\bar u_{\nu, 3/2}F^{\mu\nu}
    \f]
  */

  /*!
    \class D_DP
    \brief For decays \f$\frac12\to \frac{1}{2'}P\f$

    \f[
    {\cal M} = \bar u_{1'/2}(a+b\gamma_5)u_{1/2}
    \f]
  */
  /*!
    \class D_DV
    \brief For decays \f$\frac12\to \frac{1}{2'}V\f$

    \f[
    {\cal M} = 
    \epsilon^*_\mu\bar u_{1'/2}\left[\gamma^\mu(a+b\gamma_5)+
                                     p^\mu_{1/2}(a'+b'\gamma_5)\right]u_{1/2}
    \f]
  */
  /*!
    \class D_RP
    \brief For decays \f$\frac12\to \frac32P\f$

    \f[
    {\cal M} = \bar u^\mu_{3/2}(a+b\gamma_5)u_{1/2}p_{\mu,1/2}
    \f]
  */
  /*!
    \class D_RV
    \brief For decays \f$\frac12\to \frac32V\f$

    \f[
    {\cal M} = 
    \bar u^\mu_{3/2}\left[g_{\mu\nu}(a+b\gamma_5)+
                          p_{\mu,1/2}p_{\nu,3/2}(a'+b'\gamma_5)\right]
    u^\mu_{1/2}\epsilon^{*\nu}    
    \f]
  */
  /*!
    \class R_DP
    \brief For decays \f$\frac32\to \frac12P\f$

    \f[
    {\cal M} = \bar u_{1/2}(a+b\gamma_5)u^\mu_{3/2}p_{\mu,1/2}
    \f]
  */
  /*!
    \class R_DV
    \brief For decays \f$\frac32\to \frac12V\f$

    \f[
    {\cal M} = 
    \bar u_{1/2}\left[g_{\mu\nu}(a+b\gamma_5)+
                      p_{\mu,1/2}p_{\nu,3/2}(a'+b'\gamma_5)\right]
    u^\mu_{3/2}\epsilon^{*\nu}    
    \f]
  */
  /*!
    \class R_RP
    \brief For decays \f$\frac32\to \frac{3}{2'}P\f$

    \f[
    {\cal M} = 
    \bar u^\mu_{3'/2}\left[g_{\mu\nu}(a+b\gamma_5)+
                           p_\mu p'_\nu(a'+b'\gamma_5)\right]u^\nu_{3/2}
    \f]
  */


  /*! 
    \file Baryon_SimpleDecay_MEs.H
    For non-leptonic Hyperon decays, see hep-ph/9902351, 
    for the Omegas, cf.\ hep-ph/9905398, for \f$\Omega->\Xi^*\pi\f$, see hep-ph/0405162.
  */
};

#endif
