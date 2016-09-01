#ifndef DIRE__Tools__Splitting_H
#define DIRE__Tools__Splitting_H

#include "DIRE/Tools/Color.H"
#include "DIRE/Tools/Parton.H"
#include "DIRE/Tools/Weight.H"

namespace DIRE {

  class Kernel;

  struct Splitting {

    Parton *p_c, *p_s, *p_n, *p_l;
    Parton_Vector m_ss;
    const Kernel *p_sk;
    double m_t, m_z, m_phi, m_s, m_z2, m_phi2;
    double m_Q2, m_x, m_y;
    double m_mij2, m_mi2, m_mj2, m_mk2, m_ml2;
    double m_eta, m_t0, m_t1;
    Color_Vector m_ci, m_cj;
    int m_h[3];
    int m_type, m_cm, m_kin, m_cpl, m_kfac, m_clu;
    Weight m_w;

    Splitting(Parton *const c=NULL,Parton *const s=NULL,
	      const double &t=0.0,const double &z=0.0,
	      const double &phi=0.0):
      p_c(c), p_s(s), p_n(NULL), p_l(NULL), p_sk(NULL), 
      m_t(t), m_z(z), m_phi(phi), m_s(0.0), m_z2(0.0), m_phi2(0.0),
      m_Q2(0.0), m_x(0.0), m_y(0.0),
      m_mij2(0.0), m_mi2(0.0), m_mj2(0.0), m_mk2(0.0), m_ml2(0.0),
      m_eta(0.0), m_t0(0.0), m_t1(0.0),
      m_type(-1), m_cm(-1), m_kin(1), m_cpl(0), m_kfac(0), m_clu(0),
      m_w(0,0,0)
    { m_h[2]=m_h[1]=m_h[0]=0; }

    void SetType();

    inline void ResetCol() { m_ci.clear(); m_cj.clear(); }

    inline void AddCol(const Color &ci,const Color &cj)
    { m_ci.push_back(ci); m_cj.push_back(cj); }

  };// end of struct Splitting

  std::ostream &operator<<(std::ostream &s,const Splitting &p);

}// end of namespace DIRE

#endif
