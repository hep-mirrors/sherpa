#ifndef PHASIC__Scales__MINLO_Scale_Setter_H
#define PHASIC__Scales__MINLO_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "PDF/Main/ISR_Handler.H"

namespace PHASIC {

  struct MCS_Params {
    size_t m_i, m_j, m_k;
    ATOOLS::Flavour m_fl;
    double m_kt2, m_op2;
    ATOOLS::Vec4D m_pijt, m_pkt;
    ATOOLS::Poincare_Sequence m_lam;
    MCS_Params(const size_t &i,const size_t &j,const size_t &k,
	      const ATOOLS::Flavour &fl=ATOOLS::Flavour(kf_gluon)):
      m_i(i),m_j(j), m_k(k), m_fl(fl), m_kt2(-1.0), m_op2(-1.0) {}
    bool operator<(const MCS_Params &ck) const
    { 
      if (m_i<ck.m_i) return true; if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true; if (m_j>ck.m_j) return false;
      if (m_k<ck.m_k) return true; if (m_k>ck.m_k) return false;
      return m_fl<ck.m_fl;
    }
    void SetParams(const double &kt2,
		   const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt,
		   const ATOOLS::Poincare_Sequence &lam=
		   ATOOLS::Poincare_Sequence())
    { m_kt2=kt2; m_op2=1.0/m_kt2; m_pijt=pijt; m_pkt=pkt; m_lam=lam; }
  };// end of struct MCS_Params

  class MINLO_Scale_Setter: public Scale_Setter_Base {
  private:

    Core_Scale_Setter *p_core;
    PDF::ISR_Handler  *p_isr;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

    ATOOLS::Flavour_Vector m_f;

    ATOOLS::Cluster_Amplitude *p_ampl;

    double m_rsf, m_dr, m_q02[2], m_muf2min;
    int    m_noutmin, m_rproc, m_vproc, m_nproc, m_nfgsplit;
    int    m_cmode, m_hqmode, m_order, m_orderrs;
    int    m_usecomb, m_usepdfinfo, m_nlocpl, m_mufmode;

    PDF::CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const;

    void KT2(const ATOOLS::Cluster_Leg *li,
	     const ATOOLS::Cluster_Leg *lj,
	     const ATOOLS::Cluster_Leg *lk,MCS_Params &cs) const;

    bool Combine(ATOOLS::Cluster_Amplitude &ampl,const MCS_Params &cs) const;

    double SetScales(ATOOLS::Cluster_Amplitude *ampl,const size_t &mode);

  public:

    MINLO_Scale_Setter(const Scale_Setter_Arguments &args);

    ~MINLO_Scale_Setter();

    double Calculate(const ATOOLS::Vec4D_Vector &p,const size_t &mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

    inline ATOOLS::Cluster_Amplitude *Ampl() const { return p_ampl; }

    inline double Q02(const int i) const { return m_q02[i]; }

  };// end of class MINLO_Scale_Setter

}// end of namespace PHASIC

#endif
