#ifndef EXTRAXS__Two2Two__NSCS_Tools_H
#define EXTRAXS__Two2Two__NSCS_Tools_H

#include "ATOOLS/Math/Vector.H"

namespace NSCS {

  struct NSCS_Args {
    double m_Q, m_x1, m_x2, m_x3, m_x4, m_lam, m_phi, m_w;
    double m_e13, m_e14, m_e34, m_s14s13, m_s34s13;
    int m_sec, m_sphi;
    ATOOLS::Vec4D m_p[4], m_b, m_a;
    inline NSCS_Args(const double &Q=0.,
		     const double &x1=0.,const double &x2=0.,
		     const double &x3=0.,const double &x4=0.,
		     const double &lam=0.,const double &phi=0.):
      m_Q(Q), m_x1(x1), m_x2(x2), m_x3(x3), m_x4(x4),
      m_lam(lam), m_phi(phi), m_e13(0.), m_e14(0.), m_e34(0.),
      m_s14s13(0.), m_s34s13(0.), m_w(1.0),
      m_sec(0), m_sphi(0) { }
    inline NSCS_Args(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2,
		     const ATOOLS::Vec4D &p3,const ATOOLS::Vec4D &p4):
      m_Q(0.), m_x1(0.), m_x2(0.), m_x3(0.), m_x4(0.),
      m_lam(0.), m_phi(0.), m_e13(0.), m_e14(0.), m_e34(0.),
      m_s14s13(0.), m_s34s13(0.), m_w(0.),
      m_sec(0), m_sphi(0)
    { m_p[0]=p1; m_p[1]=p2; m_p[2]=p3; m_p[3]=p4; }
  };//end of struct NSCS_Args

  double D(const double &e3,const double &e4,const double &lam);
  double ID(const double &e3,const double &e4,const double &eta);

  double N(const double &x3,const double &x4,const double &lam);
  double IN(const double &x3,const double &x4,const double &eta);

  int ClusterNNLO(NSCS_Args &a);
  int ConstructNNLO(NSCS_Args &a,const int mode=0);

  double PSWeight(const NSCS_Args &a,const int mode=3);

  double W1314(const NSCS_Args &a,const int &mode=0);
  double W1324(const NSCS_Args &a,const int &mode=0);

  int Theta(const NSCS_Args &a,const int mode);

}// end of namespace NCSC

#endif
