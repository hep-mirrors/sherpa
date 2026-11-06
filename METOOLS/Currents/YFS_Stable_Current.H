#ifndef METOOLS__Currents__YFS_Stable_Current_H
#define METOOLS__Currents__YFS_Stable_Current_H

#include "METOOLS/Explicit/Current.H"
#include "ATOOLS/Math/DDouble.H"
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
    ATOOLS::DDouble Virtuality(const size_t &cid) const
    {
      const size_t set(cid&Full()), cnt(Count(set));
      if (m_n-cnt<cnt && Balanced(ATOOLS::Vec4D())) return SetVirtuality(Full()^set);
      return SetVirtuality(set);
    }
  };// end of struct Invariant_Table

  /*!
    Dot product of two momenta whose virtualities are known exactly.
    Vec4::LCDot forms p- = p[0]-p[3] by subtraction. 
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
    // recover whichever component cancelled, from the one that did not
    if (abs(pp)>abs(pm)) { if (pp!=DD(0.0)) pm=(p2+ptt)/pp; }
    else                 { if (pm!=DD(0.0)) pp=(p2+ptt)/pm; }
    if (abs(qp)>abs(qm)) { if (qp!=DD(0.0)) qm=(q2+qtt)/qp; }
    else                 { if (qm!=DD(0.0)) qp=(q2+qtt)/qm; }
    return DD(0.5)*(pp*qm+pm*qp)-p[1]*q[1]-p[2]*q[2];
  }

  /*!
    The YFS part of a "Y" current: its double-double momentum and
    virtuality,
  */
  class YFS_Stable_Current {
  public:
    typedef ATOOLS::DDouble DD;
  protected:
    const Invariant_Table *p_inv=nullptr;
    // double-double shadow of the current's m_p, and its virtuality
    ATOOLS::Vec4<DD> m_ph;
    DD m_p2;
    // externally supplied wide momentum (Amplitude's joint on-shell and
    // conserving projection), wins over the per-leg reconstruction
    ATOOLS::Vec4<DD> m_phx;
    bool m_hasphx=false;
    // the YFS parts of the sub-currents of the first vertex, cached
    std::vector<const YFS_Stable_Current*> m_subs;
    bool m_subsset=false;
    // the pole at a fixed scale instead of p^2, <= 0 for none
    double m_p2prop=-1.0;
    // the pole at (p + d)^2, p this current's all-outgoing momentum
    ATOOLS::Vec4D m_pshift;
    bool m_haveshift=false;
    /*
      ACRAIC: RESONANCE_EXPONENTIATION. With a non-zero coefficient a the
      propagator of this (resonant) current is multiplied by
      exp(a ln((M~^2 - q^2)/M~^2)), M~^2 = M^2 - i M Gamma, q^2 the momentum
      its pole sits at, or by the linear 1 + a ln(...) (the O(alpha) part).
    */
    double m_resexp=0.;
    bool m_reslinear=false;
    // YFS: IFI false. A vetoed line's propagator is zero.
    bool m_propveto=false;

    //! the wide momentum of an external current of momentum p
    void SetPWide(const ATOOLS::Vec4D &p)
    {
      if (m_hasphx) { m_ph=m_phx; return; }
      DD pp(DD(p[0])+DD(p[3])), pm(DD(p[0])-DD(p[3]));
      const DD pt2(DD(p[1])*DD(p[1])+DD(p[2])*DD(p[2]));
      if (abs(pp)>abs(pm)) { if (pp!=DD(0.0)) pm=(m_p2+pt2)/pp; }
      else                 { if (pm!=DD(0.0)) pp=(m_p2+pt2)/pm; }
      m_ph=ATOOLS::Vec4<DD>(DD(0.5)*(pp+pm),DD(p[1]),DD(p[2]),DD(0.5)*(pp-pm));
    }
    //! an external current of momentum p and mass m
    void SetExternal(const ATOOLS::Vec4D &p,const double &m)
    { m_p2=DD(m*m); SetPWide(p); }

    //! Current::SumMomenta: p from the sub-currents js, its virtuality as
    //! sum p_i^2 + 2 sum p_i.p_j, narrowed to double once, at the end
    void SumWideMomenta(const Current_Vector &js,ATOOLS::Vec4D &p)
    {
      if (!m_subsset) {
        m_subs.resize(js.size());
        for (size_t i(0);i<js.size();++i)
          m_subs[i]=dynamic_cast<const YFS_Stable_Current*>(js[i]);
        m_subsset=true;
      }
      std::vector<ATOOLS::Vec4<DD> > ph(js.size());
      std::vector<DD> p2(js.size());
      for (size_t i(0);i<js.size();++i) {
        if (m_subs[i]) { ph[i]=m_subs[i]->m_ph; p2[i]=m_subs[i]->m_p2; }
        else {
          // a current without YFS part: its double momentum
          const ATOOLS::Vec4D q(js[i]->P());
          ph[i]=ATOOLS::Vec4<DD>(DD(q[0]),DD(q[1]),DD(q[2]),DD(q[3]));
          p2[i]=DD(q.Abs2());
        }
      }
      m_ph=ATOOLS::Vec4<DD>();
      m_p2=DD(0.0);
      for (size_t i(0);i<js.size();++i) {
        m_p2+=p2[i];
        for (size_t j(i+1);j<js.size();++j)
          m_p2+=DD(2.0)*YFSExactLCDot(ph[i],p2[i],ph[j],p2[j]);
        m_ph+=ph[i];
      }
      p=ATOOLS::Vec4D((double)m_ph[0],(double)m_ph[1],
                      (double)m_ph[2],(double)m_ph[3]);
    }

    /*!
      YFS: COMIX_STABLE_PROPAGATORS. The virtuality p2 and momentum ph the
      propagator of current c uses, shifted: whether the line carries an
      ACRAIC shift (then ph is the shifted momentum). Returns false without
      a valid invariant table.
    */
    bool StableLine(const Current &c,DD &p2,ATOOLS::Vec4<DD> &ph,
                    bool &shifted) const
    {
      shifted=m_haveshift;
      if (p_inv==nullptr || !p_inv->m_valid) return false;
      if (!shifted) { p2=p_inv->Virtuality(c.CId()); ph=m_ph; return true; }
      const ATOOLS::Vec4D &dv(m_pshift);
      if (p_inv->Balanced(dv)) {
        const size_t sbar(p_inv->Full()^(c.CId()&p_inv->Full()));
        const ATOOLS::Vec4<DD> q(p_inv->SetMomentum(sbar));
        ph=ATOOLS::Vec4<DD>(-q[0],-q[1],-q[2],-q[3]);
        p2=p_inv->SetVirtuality(sbar);
        return true;
      }
      const ATOOLS::Vec4<DD> d{DD{dv[0]},DD{dv[1]},DD{dv[2]},DD{dv[3]}};
      const DD d2(d[0]*d[0]-d[1]*d[1]-d[2]*d[2]-d[3]*d[3]);
      const ATOOLS::Vec4<DD> &h(m_ph);
      ph=ATOOLS::Vec4<DD>(h[0]+d[0],h[1]+d[1],h[2]+d[2],h[3]+d[3]);
      p2=m_p2+DD(2.0)*YFSExactLCDot(h,m_p2,d,d2)+d2;
      return true;
    }

  public:
    YFS_Stable_Current(const double &m): m_p2(m*m) {}
    virtual ~YFS_Stable_Current() {}

    inline void SetInvariantTable(const Invariant_Table *const inv) { p_inv=inv; }
    inline void SetPWideExt(const ATOOLS::Vec4<DD> &ph)
    { m_phx=ph; m_hasphx=true; }

    inline void SetPropScale(const double &s) { m_p2prop=s; }
    inline void SetPropShift(const ATOOLS::Vec4D &d)
    { m_pshift=d; m_haveshift=(d[0]!=0.||d[1]!=0.||d[2]!=0.||d[3]!=0.); }
    inline void SetResonanceExponent(const double a,const bool linear=false)
    { m_resexp=a; m_reslinear=linear; }
    inline void SetPropVeto(const bool v) { m_propveto=v; }

    inline double PropScale() const { return m_p2prop; }
    inline const ATOOLS::Vec4D &PropShift() const { return m_pshift; }
    inline bool HavePropShift() const { return m_haveshift; }
    inline double ResonanceExponent() const { return m_resexp; }
    inline bool PropVeto() const { return m_propveto; }
  };// end of class YFS_Stable_Current

}// end of namespace METOOLS

#endif
