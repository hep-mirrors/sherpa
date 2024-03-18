#ifndef ALARIC__Tools__Splitting_H
#define ALARIC__Tools__Splitting_H

#include "ALARIC++/Tools/Color.H"
#include "ALARIC++/Tools/Parton.H"
#include "ALARIC++/Tools/Weight.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Poincare.H"

namespace ALARIC {

  class Kernel;

  struct Splitting {

    Parton *p_c, *p_s, *p_n, *p_l;
    const Kernel *p_sk;
    double m_t, m_z, m_phi, m_s, m_z2, m_phi2;
    double m_q2, m_Q2, m_x, m_y, m_kt2, m_kap;
    double m_mij2, m_mi2, m_mj2, m_mk2, m_ml2;
    double m_eta, m_t0, m_t1, m_J;
    Color_Vector m_ci, m_cj;
    ATOOLS::Vec4D m_pi, m_pj, m_pl, m_pk, m_Kt, m_K;
    ATOOLS::Vec4D_Vector m_p;
    ATOOLS::Poincare_Sequence m_lam;
    int m_h[3], m_mode, m_soft, m_mec, m_iink;
    int m_type, m_cm, m_cpl, m_kfac, m_clu;
    Weight m_w;
    std::vector<double> m_vars;
    std::vector<int> m_rcl;

    Splitting(Parton *const c=NULL,Parton *const s=NULL,
	      const double &t=0.0,const double &z=0.0,
	      const double &phi=0.0):
      p_c(c), p_s(s), p_n(NULL), p_l(NULL), p_sk(NULL), 
      m_t(t), m_z(z), m_phi(phi), m_s(0.0), m_z2(0.0), m_phi2(0.0),
      m_Q2(0.0), m_x(0.0), m_y(0.0), m_kt2(0.0),
      m_mij2(0.0), m_mi2(0.0), m_mj2(0.0), m_mk2(0.0), m_ml2(0.0),
      m_eta(0.0), m_t0(0.0), m_t1(0.0), m_J(1.0), m_mode(0), m_soft(0),
      m_mec(0), m_iink(0), m_type(-1), m_cm(-1), m_cpl(0), m_kfac(0),
      m_clu(0), m_w(0,0,0)
    { m_h[2]=m_h[1]=m_h[0]=0; }

    void SetType();

    inline void ResetCol() { m_ci.clear(); m_cj.clear(); }

    inline void AddCol(const Color &ci,const Color &cj)
    { m_ci.push_back(ci); m_cj.push_back(cj); }

  };// end of struct Splitting

  std::ostream &operator<<(std::ostream &s,const Splitting &p);

}// end of namespace ALARIC

#endif
