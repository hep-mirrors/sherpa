
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"

#include "ATOOLS/Org/Message.H"

namespace HIGGS {

  extern ATOOLS::Function_Base *s_mt, *s_mb, *s_mc;

  extern double alpha0, m_W, m_Z, s2W, c2W, sW, cW;
  extern double m_u, m_d, m_s;
  extern double m_e, m_mu, m_tau;
  extern double e_u, e_c, e_t, e_d, e_s, e_b;
  extern double I_3u, I_3c, I_3t, I_3d, I_3s, I_3b;

  extern double G_F;

  const double N_c=3.0, C_A=3.0;
  const double C_F=4.0/3.0;
  const double N_f=5.0;

  extern double sumQsq, sumQ4, sumQ6;

  extern double mu_sq;

  const Complex I(0.0,1.0);

  extern AMEGIC::Basic_Sfuncs *s_bs;

  inline Complex spa(int i,int j) { return s_bs->S0(i-1,j-1); }
  inline Complex spb(int i,int j) { return s_bs->S1(i-1,j-1); }

  inline double sij(int i,int j)
  { return (spa(i,j)*spb(j,i)).real(); }
}

using namespace HIGGS;

inline double M_t(const double &mu)
{
  return (*s_mt)(mu*mu);
}

inline double M_b(const double &mu)
{
  return (*s_mb)(mu*mu);
}

inline double M_c(const double &mu)
{
  return (*s_mc)(mu*mu);
}

inline double alpha_s(const double &mur)
{
  return (*MODEL::as)(mur*mur);
}

inline double alpha_QED(const double &mur)
{
  return (*MODEL::aqed)(mur*mur);
}
