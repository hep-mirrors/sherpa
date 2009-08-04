#ifndef HADRONS_ME_Library_Eta_Decay_MEs_H
#define HADRONS_ME_Library_Eta_Decay_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {
  class Eta_PVV : public HD_ME_Base {
    int     m_npol1,m_npol2,m_ff;
    double  m_mrho2,m_momega2,m_mrho,m_momega,m_Grho,m_Gomega,m_fP,m_g_rhoG;
    Complex m_global;
    Complex D(const ATOOLS::Vec4D *);
    Complex E(const ATOOLS::Vec4D *);
  public:
    Eta_PVV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *);
    void SetModelParameters(GeneralModel);
  };

  class Eta_PPV : public HD_ME_Base {
    int     m_npol,m_ff;
    double  m_VDM_mass,m_VDM_width,m_fP;
    double  m_mpipi2,m_mrho2,m_pref;
    Complex m_global;
    Complex Formfactor(const double);
    Complex Omnes_Formfactor(const double,const double);
  public:
    Eta_PPV(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name), 
      m_npol(3), m_ff(0), m_VDM_mass(0.), m_VDM_width(0.), m_fP(0.) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *);
    void SetModelParameters(GeneralModel);
  };

  class Eta_PPP : public HD_ME_Base {
  public:
    Eta_PPP(ATOOLS::Flavour * flavs,int n,int* indices,std::string name) :
      HD_ME_Base(flavs,n,indices,name) {};
    void operator()(const ATOOLS::Vec4D *,METOOLS::Spin_Amplitudes *);
    void SetModelParameters(GeneralModel);
  };

  /*!
    \class Eta_PVV
    \brief For decays \f$\eta\to\pi\gamma\gamma\f$

    \f[
    \f]
    cf.\ hep-ph/0112150.
   */
  
  /*!
    \class Eta_PPV
    \brief For decays \f$\eta\to\pi\pi\gamma\f$

    \f[
    \f]
    cf.\ hep-ph/0112150.
   */
  
  /*!
    \class Eta_PPP
    \brief For decays \f$\eta\to\pi\pi\pi\f$ and \f$\eta'\to\eta\pi\pi\f$

    \f[
    \f]
    cf.\ hep-ph/0301058 and hep-ph/0511076.
   */
  

}

#endif
