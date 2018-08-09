#ifndef CFPSHOWER_Shower_Splitting_Parameters_H
#define CFPSHOWER_Shower_Splitting_Parameters_H

//#include "CFPSHOWER/Tools/Splitting.H"

namespace CFPSHOWER {
  class Kernel;

  struct Splitting_Parameters {
    double   m_kt2, m_z, m_y, m_phi, m_x, m_wk, m_ws, m_q2, m_mu2;
    int      m_mode, m_kin, m_col, m_kmode;
    Kernel * p_sf;
    size_t   m_idi, m_idj, m_idk;
    //ATOOLS::Vec4D m_pbt;
    CS_Parameters(const double & kt2, const double & z,
		  const double & y,   const double & phi,
		  const double & x,   const double & q2,
		  const int mode=-1,  const int kin=0,  const int kmode=0):
      m_kt2(kt2), m_z(z), m_y(y), m_phi(phi), m_x(x), m_wk(1.0), m_ws(kt2),
      m_q2(q2), m_mu2(kt2), m_mode(mode), m_kin(kin), m_col(0), m_kmode(kmode), p_sf(NULL),
      m_idi(0), m_idj(0), m_idk(0) {}
  };// end of struct Splitting_Parameters

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs);

}// end of namespace CFPSHOWER

#endif
