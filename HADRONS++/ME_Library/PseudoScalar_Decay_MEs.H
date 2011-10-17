#ifndef HADRONS_ME_Library_PseudoScalar_Decay_MEs_H
#define HADRONS_ME_Library_PseudoScalar_Decay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Math/MyComplex.H"

namespace HADRONS {
  class P_PP : public HD_ME_Base {
    Complex m_global;
  public:
    P_PP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_LNu : public HD_ME_Base {
    Complex m_global;
  public:
    P_LNu(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_PV : public HD_ME_Base {
    int     m_npol;
    Complex m_global;
  public:
    P_PV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_npol(3), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_VV : public HD_ME_Base {
    int     m_npol1,m_npol2;
    Complex m_global;
  public:
    P_VV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_npol1(3), m_npol2(3), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class S_VV : public HD_ME_Base {
    int     m_npol1,m_npol2;
    Complex m_global;
  public:
    S_VV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_npol1(3), m_npol2(3), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_PT : public HD_ME_Base {
    Complex m_global;
  public:
    P_PT(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_PPP : public HD_ME_Base {
    int     m_ff;
    double  m_g,m_h,m_j,m_k,m_f;
    Complex m_global;
    double  Formfactor(const ATOOLS::Vec4D *);
  public:
    P_PPP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), 
      m_ff(0), m_g(0.), m_h(0.), m_j(0.), m_k(0.), m_f(0.),
      m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_VFF : public HD_ME_Base {
    int     m_npol;
    bool    m_VDM;
    double  m_VDM_mass, m_VDM_width;
    Complex m_global;
  public:
    P_VFF(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), 
      m_npol(3), m_VDM(false), m_VDM_mass(0.), m_VDM_width(0.), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_PLNu : public HD_ME_Base {
    int     m_ff;
    double  m_Norm,m_fP,m_aplus0,m_aplus1,m_azero0,m_azero1,m_masssqr_diff;
    Complex m_global;

    double fplus(const double,const double,const double);
    double fzero(const double,const double,const double);
  public:
    P_PLNu(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_ff(0),m_fP(1.),m_Norm(1.),
      m_aplus0(0.),m_aplus1(0.),m_azero0(0.),m_azero1(0.),m_masssqr_diff(0.),
      m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_PPLNu : public HD_ME_Base {
    int     m_ff;
    double  m_fP,m_f1_0,m_lambda1,m_f2_0,m_lambda2,m_g_0,m_kappa;
    Complex m_global;

    double f1(const double,const double,const double);
    double f2(const double,const double,const double);
    double g(const double,const double,const double);
  public:
    P_PPLNu(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_ff(0),m_fP(1.),
      m_f1_0(1.),m_lambda1(0.),m_f2_0(1.),m_lambda2(0.),m_g_0(1.),m_kappa(0.),
      m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  class P_FFFF : public HD_ME_Base {
    bool m_exchange;
    Complex m_global;
  public:
    P_FFFF(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), m_exchange(false), m_global(Complex(0.,0.)) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters(GeneralModel);
  };

  /*!
    \class P_PP
    \brief For decays \f$P(S)\to P_1P_2\f$, \f$P(S)\to S_1S_2\f$

    \f[
    {\cal M} = g\\exp(i\\phi)
    \f]
  */
  /*!
    \class P_LNu
    \brief For decays \f$P\to \ell\bar\nu_l\f$

    \f[
    {\cal M} = \frac{G_F}{\sqrt{2}} f_Pm_P V_{\rm CKM}
    \bar u_\ell(1-\gamma_5)u_{\bar\nu}
    \f]
  */
  /*!
    \class P_PV
    \brief For decays \f$P\to P'V\f$ and \f$S\to S'V\f$.

    \f[
    {\cal M} = g\epsilon^\mu_V(p_P+p_{P'})_\mu
    \f]
  */
  /*!
    \class P_VV
    \brief For decays \f$P\to VV'\f$

    \f[
    {\cal M} = g\epsilon_{\mu\nu\rho\sigma}
    p^\mu_V\epsilon^\nu_V p^\rho_{V'}\epsilon^\sigma_{V'}
    \f]
  */
  /*!
    \class S_VV
    \brief For decays \f$S\to VV'\f$

    \f[
    {\cal M} = g\epsilon^\mu_V\epsilon_{\mu, V'}
    \f]
  */
  /*!
    \class P_PT
    \brief For decays \f$P\to P'T\f$

    \f[
    {\cal M} = g\epsilon^{\mu\nu}_Tp_{\mu, P}p_{\nu, P'}
    \f]
  */
  /*!
    \class P_VFF
    \brief For decays \f$P\to V\ell\bar\ell\f$

    \f[
    {\cal M} = \frac{g}{q^2_{\ell\bar\ell}}
    \epsilon_{\mu\nu\rho\sigma}p^\mu_V\epsilon^\nu_Vq^\rho_{\ell\bar\ell}
    (\bar u_\ell\gamma^\sigma u_{\bar\ell})
    \f]

    In case we want to do a VDM model, the propagator term \f$1/q^2\f$ should be replaced 
    with \f$-(m_{VDM}^2-i\Gamma_{VDM}m_{VDM})/(q^2-m_{VDM}^2+i\Gamma_{VDM}m_{VDM})\f$.
  */
  /*!
    \class P_PLNu
    \brief For decays \f$P\to P'\ell\bar\nu_l\f$

    \f[
    {\cal M} = \frac{G_F}{\sqrt{2}}N_{PP'}
    \left(f_+(q_{\ell\nu}^2)(p+p')_\mu+f_-(q_{\ell\nu}^2)(p-p')_\mu\right)
    \bar u_\ell\gamma^\mu(1-\gamma_5)u_{\bar\nu}
    \f]\,,
    where the form factors are given through the expansion
    \f[
    f_+(p^2,{p'}^2,q^2) &=& a_+^0 + a_+^1\frac{q^2}{f_{P'}^2}\nonumber
    f_-(p^2,{p'}^2,q^2) &=& a_-^0 + a_-^1\frac{p^2-{p'}^2}{f_{P'}^2}\,.
    \f]
    By default, these values are set to \f$a_+^0 = -1\f$ and $a_-^0=a_-^1=a_+^1=0\f$.
  */
  /*!
    \class P_FFFF
    \brief For decays \f$P\to \ell\bar\ell\ell'\bar\ell'\f$

    \f[
    {\cal M} = \frac{g}{q^2_{\ell\bar\ell}q^2_{\ell'\bar\ell'}}
    \epsilon_{\mu\nu\rho\sigma}
    (\bar u_\ell\gamma^\mu u_{\bar\ell})q^\nu_{\ell\bar\ell}
    (\bar u_\ell'\gamma^\rho u_{\bar\ell'})q^\sigma_{\ell'\bar\ell'}
    \f]

    In case we want to do a VDM model, the propagator terms \f$1/q^2\f$ should
    be replaced with \f$-(m_{VDM}^2-i\Gamma_{VDM}m_{VDM})/
    (q^2-m_{VDM}^2+i\Gamma_{VDM}m_{VDM})\f$.  Also, exchange diagrams for 
    four identical leptons should be considered.
  */
  

};

#endif
