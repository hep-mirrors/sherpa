#ifndef METOOLS__Currents__YFS_Stable_Current_H
#define METOOLS__Currents__YFS_Stable_Current_H

/*
  YFS: COMIX_STABLE_PROPAGATORS (arXiv:2406.07671, sec. II A).

  Comix currents for the real-emission processes YFS builds (YFS_Process
  marks them with "YFSStable" in Process_Info::m_special; COMIX::Amplitude
  then asks for the "Y" currents of this header's classes instead of "D").
  Only the propagators differ, and only in where their p^2 comes from:

  - the external invariants s_ab = 2 p_a.p_b are evaluated without
    cancellation (Invariant_Table, filled by COMIX::Amplitude);
  - a line whose ACRAIC propagator shift d puts back the photons the
    partition legs omit (p_S + d = -p_Sbar) takes p^2 and its momentum from
    the complementary legs: for an emission line, m^2 + 2 p.k exactly;
  - otherwise a shifted line is evaluated in double-double, and an unshifted
    current uses the complementary legs only when the legs balance.

  Every other Comix process uses the base classes unchanged.
*/

#include "METOOLS/Explicit/Current.H"
#include <cmath>

namespace METOOLS {

  struct Invariant_Table {
    bool m_valid=false;
    size_t m_n=0;
    std::vector<ATOOLS::DDouble> m_m2, m_s;   // m_s[a*m_n+b], a<b
    std::vector<ATOOLS::Vec4D> m_p;           // the legs, all outgoing
    ATOOLS::Vec4D m_R;                        // sum of all legs
    double m_tol=0.;                          // balance tolerance [GeV]
    inline size_t Full() const { return (size_t(1)<<m_n)-1; }
    inline size_t Count(const size_t &set) const
    { size_t c(0); for (size_t a(0);a<m_n;++a) if (set&(size_t(1)<<a)) ++c; return c; }
    // R + d = 0 within m_tol: the legs plus the extra momentum d balance
    inline bool Balanced(const ATOOLS::Vec4D &d) const
    {
      for (int k(0);k<4;++k) if (std::abs(m_R[k]+d[k])>m_tol) return false;
      return true;
    }
    // p^2 of the leg set with bit mask set, from the pair invariants
    ATOOLS::DDouble SetVirtuality(const size_t &set) const
    {
      ATOOLS::DDouble p2(0.0);
      for (size_t a(0);a<m_n;++a) {
        if (!(set&(size_t(1)<<a))) continue;
        p2+=m_m2[a];
        for (size_t b(a+1);b<m_n;++b)
          if (set&(size_t(1)<<b)) p2+=m_s[a*m_n+b];
      }
      return p2;
    }
    // the momentum of the leg set, in double-double
    ATOOLS::Vec4<ATOOLS::DDouble> SetMomentum(const size_t &set) const
    {
      typedef ATOOLS::DDouble DD;
      DD q[4]={DD(0.0),DD(0.0),DD(0.0),DD(0.0)};
      for (size_t a(0);a<m_n;++a)
        if (set&(size_t(1)<<a)) for (int k(0);k<4;++k) q[k]+=DD(m_p[a][k]);
      return ATOOLS::Vec4<DD>(q[0],q[1],q[2],q[3]);
    }
    // p^2 of an UNSHIFTED current: the complement only if the legs balance
    // and it is the smaller set
    ATOOLS::DDouble Virtuality(const size_t &cid) const
    {
      const size_t set(cid&Full()), cnt(Count(set));
      if (m_n-cnt<cnt && Balanced(ATOOLS::Vec4D())) return SetVirtuality(Full()^set);
      return SetVirtuality(set);
    }
  };// end of struct Invariant_Table

  /*!
    Dot product of two momenta whose virtualities are known exactly: the
    light-cone component that cancels is recovered from p+ p- = p^2 + pT^2
    (as Current.C's ExactLCDot).
  */
  inline ATOOLS::DDouble YFSExactLCDot(const ATOOLS::Vec4<ATOOLS::DDouble> &p,
                                       const ATOOLS::DDouble &p2,
                                       const ATOOLS::Vec4<ATOOLS::DDouble> &q,
                                       const ATOOLS::DDouble &q2)
  {
    typedef ATOOLS::DDouble DD;
    const DD ptt(p[1]*p[1]+p[2]*p[2]), qtt(q[1]*q[1]+q[2]*q[2]);
    DD pp(p[0]+p[3]), pm(p[0]-p[3]);
    DD qp(q[0]+q[3]), qm(q[0]-q[3]);
    if (abs(pp)>abs(pm)) { if (pp!=DD(0.0)) pm=(p2+ptt)/pp; }
    else                 { if (pm!=DD(0.0)) pp=(p2+ptt)/pm; }
    if (abs(qp)>abs(qm)) { if (qp!=DD(0.0)) qm=(q2+qtt)/qp; }
    else                 { if (qm!=DD(0.0)) qp=(q2+qtt)/qm; }
    return DD(0.5)*(pp*qm+pm*qp)-p[1]*q[1]-p[2]*q[2];
  }

  /*!
    The YFS line: the virtuality p2 and momentum ph its propagator uses.
    shifted: whether the line carries an ACRAIC shift (then ph is the
    shifted momentum). Returns false without a valid table.
  */
  class YFS_Stable_Current {
  protected:
    const Invariant_Table *p_inv=nullptr;
  public:
    virtual ~YFS_Stable_Current() {}
    inline void SetInvariantTable(const Invariant_Table *const inv) { p_inv=inv; }
    bool StableLine(const Current &c,ATOOLS::DDouble &p2,
                    ATOOLS::Vec4<ATOOLS::DDouble> &ph,bool &shifted) const
    {
      typedef ATOOLS::DDouble DD;
      shifted=c.HavePropShift();
      if (p_inv==nullptr || !p_inv->m_valid) return false;
      if (!shifted) { p2=p_inv->Virtuality(c.CId()); ph=c.PH(); return true; }
      const ATOOLS::Vec4D &dv(c.PropShift());
      if (p_inv->Balanced(dv)) {
        const size_t sbar(p_inv->Full()^(c.CId()&p_inv->Full()));
        const ATOOLS::Vec4<DD> q(p_inv->SetMomentum(sbar));
        ph=ATOOLS::Vec4<DD>(-q[0],-q[1],-q[2],-q[3]);
        p2=p_inv->SetVirtuality(sbar);
        return true;
      }
      const ATOOLS::Vec4<DD> d{DD{dv[0]},DD{dv[1]},DD{dv[2]},DD{dv[3]}};
      const DD d2(d[0]*d[0]-d[1]*d[1]-d[2]*d[2]-d[3]*d[3]);
      const ATOOLS::Vec4<DD> &h(c.PH());
      ph=ATOOLS::Vec4<DD>(h[0]+d[0],h[1]+d[1],h[2]+d[2],h[3]+d[3]);
      p2=c.P2H()+DD(2.0)*YFSExactLCDot(h,c.P2H(),d,d2)+d2;
      return true;
    }
  };// end of class YFS_Stable_Current

}// end of namespace METOOLS

#endif
