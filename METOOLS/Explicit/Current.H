#ifndef METOOLS_Explicit_Current_H
#define METOOLS_Explicit_Current_H

#include "METOOLS/Explicit/C_Object.H"
#include "METOOLS/Main/Polarization_Index.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/DDouble.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Node.H"

#include <vector>

namespace METOOLS {

  class Vertex;

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  typedef std::vector<size_t> SizeT_Vector;

  typedef std::vector<std::string>   String_Vector;
  typedef std::vector<String_Vector> String_Matrix;

  typedef ATOOLS::Node<std::string> Graph_Node;

  typedef std::vector<Vertex*> Vertex_Vector;

  char ParticleType(const ATOOLS::Flavour &fl);

  struct Current_Key {
    ATOOLS::Flavour m_fl;
    MODEL::Model_Base *p_model;
    size_t m_n;
    inline Current_Key(const ATOOLS::Flavour &fl,
		       MODEL::Model_Base *const model,
		       const size_t &n): 
      m_fl(fl), p_model(model), m_n(n) {}
    std::string Type() const;
  };// end of struct Current_Key

  class Current {
  protected:

    ATOOLS::Flavour m_fl;

    Vertex_Vector m_in, m_out;
    Int_Vector    m_id, m_fid;

    ATOOLS::Vec4D  m_p;
    // Double-double shadow of m_p, used ONLY for the p^2 bookkeeping below.
    // m_p itself stays double on purpose: at the leaves it holds a phase-space
    // momentum that is exactly a pair of doubles, so widening it there would
    // add digits that carry no information. Above the leaves m_p is a SUM, and
    // the sum is where the digits go missing - an internal line built from
    // nearly-collinear externals reaches an O(E^2 theta^2) virtuality through
    // O(E^2) partial sums. Carrying that sum at ~106 bits costs two doubles
    // per current on a layer that is a negligible share of the flops.
    ATOOLS::Vec4<ATOOLS::DDouble> m_ph;
    // Externally supplied wide momentum (Amplitude's joint on-shell /
    // conserving projection). When set it WINS over the per-leg
    // reconstruction in SetPWide(), because only a joint projection can
    // satisfy sum(p)=0 and p^2=m^2 at the same time.
    ATOOLS::Vec4<ATOOLS::DDouble> m_phx;
    bool m_hasphx;

    ATOOLS::DDouble m_p2;
    CObject_Matrix m_j;

    Polarization_Index m_h;
    SizeT_Vector       m_hm;

    size_t m_key, m_cid, m_ntc;
    char   m_type;

    std::vector<int> m_order;

    double m_mass, m_width;
    bool   m_msv, m_zero;
    int    m_dir, m_cut, m_osd;

    Current *p_sub;

    mutable std::string m_psinfo;

    virtual std::string CLabel() const;

    void CollectGraphs(Graph_Node *graph,
		       const std::string &lastv) const;

  public:

    Current(const Current_Key &key);

    virtual ~Current();

    // member functions
    void FindPermutations();
    void InitPols(const Int_Vector &pols);

    void AddJ(CObject *const j);

    virtual void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			    const int cr,const int ca,const int mode) = 0;
    virtual void SetGauge(const ATOOLS::Vec4D &k) = 0;

    virtual void AddPropagator() = 0;

    template <class SType> inline void 
    Contract(const Current &c,const Int_Vector &pols,
	     std::vector<std::complex<SType> > &ress,
	     const size_t &offset=0) const;

    virtual std::string Format(const CObject *c) const = 0;

    virtual char Type() const = 0;

    void SetId(const Int_Vector &id);

    void CollectGraphs(Graph_Node *graph) const;

    void DetachOut(Vertex *const v);

    std::string PSInfo() const;

    void ResetJ();
    void ResetZero();
    void Evaluate();

    void Print() const;

    // inline functions
    inline const CObject_Matrix &J() const { return m_j; }

    // External / dipole-mapped legs are on their mass shell by construction,
    // so seed the virtuality with the exact m^2 rather than with p.Abs2():
    // the latter carries an O(E^2 eps) error that would swamp the small
    // denominators this bookkeeping exists to protect.
    /*!
      Rebuild the WIDE momentum m_ph exactly on the mass shell given by m_p2.

      The stored double components satisfy p^2 = m^2 only to about
      2 E ulp(E) ~ 6e-13 GeV^2. That floor is ABSOLUTE: it comes from the
      representation of an O(E) component, so it does not shrink with the mass
      and no amount of arithmetic precision recovers it.

      Everything downstream takes the VIRTUALITY from m_p2 (the exact nominal
      m^2; ExactLCDot reconstructs the cancelling light-cone component from it)
      but the KINEMATICS from m_ph. Leaving m_ph a widened copy of the raw
      double makes those two mutually inconsistent at 6e-13, and that
      inconsistency is inherited by every current built from this leg. It stays
      invisible until a propagator denominator falls below it -- the collinear
      ISR region, 2 p.k <~ 1e-12 -- which is exactly where the real ME departed
      from the eikonal theorem.

      Keep the large light-cone component and pT exactly as given and solve the
      small one from p+ p- = p^2 + pT^2, the identity ExactLCDot already uses.
      The momentum moves by ~3e-15 in one component; the shell becomes exact.
    */
    inline void SetPWideExt(const ATOOLS::Vec4<ATOOLS::DDouble> &ph)
    { m_phx=ph; m_hasphx=true; }

    inline void SetPWide()
    {
      if (m_hasphx) { m_ph=m_phx; return; }
      { extern long g_pwext_miss; if (++g_pwext_miss==1)
          std::cerr<<"@@@ PWMISS first hit, flav="<<m_fl<<" cid="<<m_cid<<std::endl; }
      typedef ATOOLS::DDouble DD;
      DD pp(DD(m_p[0])+DD(m_p[3])), pm(DD(m_p[0])-DD(m_p[3]));
      const DD pt2(DD(m_p[1])*DD(m_p[1])+DD(m_p[2])*DD(m_p[2]));
      if (abs(pp)>abs(pm)) { if (pp!=DD(0.0)) pm=(m_p2+pt2)/pp; }
      else                 { if (pm!=DD(0.0)) pp=(m_p2+pt2)/pm; }
      m_ph=ATOOLS::Vec4<DD>(DD(0.5)*(pp+pm),DD(m_p[1]),DD(m_p[2]),
                            DD(0.5)*(pp-pm));
    }

    inline void SetP(const ATOOLS::Vec4D &p)
    { m_p=p; m_p2=ATOOLS::sqr(m_mass); SetPWide(); }

    // For the few call sites that copy a momentum off another current, where
    // the on-shell assumption above need not hold: carry that current's
    // virtuality across explicitly instead of re-deriving it.
    inline void SetP(const ATOOLS::Vec4D &p,const double &p2)
    { m_p=p; m_ph=p; m_p2=p2; }

    inline void SetFId(const Int_Vector &fid) { m_fid=fid; }
    inline void SetKey(const size_t &key)     { m_key=key; }

    inline void SetOrder(const std::vector<int> &o) { m_order=o; }
    inline void SetNTChannel(const size_t &ntc) { m_ntc=ntc;   }

    inline void SetDirection(const int &dir) { m_dir=dir; }
    inline void SetCut(const int &cut)       { m_cut=cut; }
    inline void SetOnShell(const int &osd)   { m_osd=osd; }

    inline void SetSub(Current *const sub) { p_sub=sub; }

    inline ATOOLS::Flavour Flav() const { return m_fl; }

    inline ATOOLS::Flavour RFlav() const 
    { return m_fl.IsAnti()?m_fl.Bar():m_fl; }
    inline ATOOLS::Flavour AFlav() const 
    { return m_fl.IsAnti()?m_fl:m_fl.Bar(); }

    inline ATOOLS::Vec4D P() const { return m_p; }

    // Narrowed for the ~60 call sites that only need a double; the recursion
    // itself must use PH()/P2H() so the accumulation stays wide.
    inline double P2() const { return (double)m_p2; }

    inline const ATOOLS::Vec4<ATOOLS::DDouble> &PH() const { return m_ph; }

    /*!
      The on-shell wide momentum, narrowed to the scalar type a current tower
      is instantiated at. Leaf currents must build their spinors and
      polarisations from THIS, not from the raw Vec4D they were handed: the
      raw one is off its nominal mass shell by ~6e-13 GeV^2 (see SetPWide),
      and a spinor built from it is inconsistent with the m^2 passed alongside
      it. At SType=double this is the raw momentum back again to within an
      ulp, so it costs nothing there; at the wider PMODEs it is what lets the
      extra mantissa actually carry mass information instead of re-rounding a
      momentum that was never on shell.
    */
    template <typename SType>
    inline ATOOLS::Vec4<SType> PW() const
    { return ATOOLS::Vec4<SType>((SType)m_ph[0],(SType)m_ph[1],
                                 (SType)m_ph[2],(SType)m_ph[3]); }
    inline const ATOOLS::DDouble &P2H() const { return m_p2; }

    inline const Int_Vector &Id() const   { return m_id;   }
    inline const Int_Vector &FId() const  { return m_fid;  }
    inline size_t            CId() const  { return m_cid;  }
    inline size_t            Key() const  { return m_key;  }

    inline const std::vector<int> &Order() const { return m_order; }
    inline int Order(const size_t &id) const { return m_order[id]; }

    inline size_t NTChannel() const { return m_ntc;  }

    inline int Direction() const { return m_dir; }
    inline int Cut() const       { return m_cut; }
    inline int OnShell() const   { return m_osd; }

    inline Current *Sub() const { return p_sub; }

    inline void AttachIn(Vertex *const v)  { m_in.push_back(v);  }
    inline void AttachOut(Vertex *const v) { m_out.push_back(v); }

    inline size_t NIn() const  { return m_in.size();  }
    inline size_t NOut() const { return m_out.size(); }

    inline const Vertex_Vector &In() const  { return m_in;  }
    inline const Vertex_Vector &Out() const { return m_out; }

    inline bool Dangling() const { return m_out.empty(); }

    inline double Mass() const  { return m_mass;  }
    inline double Width() const { return m_width; }

    inline bool Zero() const { return m_zero; }

    inline const Polarization_Index &H() const { return m_h; }

    inline SizeT_Vector &HM() { return m_hm; }

  };// end of class Current
  
  template <typename SType>
  class Current_Contractor {
  public:

    virtual ~Current_Contractor() {}

    // member functions
    virtual void 
    SContract(const Current &c,const Int_Vector &pols,
	      std::vector<std::complex<SType> > &ress,
	      const size_t &offset=0) const = 0;

  };// end of class Current

  // Amplitude keeps its result containers in double whatever PMODE is, and
  // Current::Contract dispatches on the element type of those containers. A
  // current tower instantiated at a wider SType therefore has to answer a
  // Contract<double> call as well: this mixin does the contraction sum at
  // SType and rounds only the finished amplitude down to double, so no
  // cancellation is ever exposed to double precision. For SType == double it
  // is empty - the class already is a Current_Contractor<double>.
  template <typename Derived,typename SType>
  class Downcast_Contractor: public Current_Contractor<double> {
  public:

    void SContract(const Current &c,const Int_Vector &pols,
                   std::vector<std::complex<double> > &ress,
                   const size_t &offset=0) const
    {
      std::vector<std::complex<SType> >
        tmp(ress.size(),std::complex<SType>(0.0,0.0));
      static_cast<const Derived*>(this)->SContract(c,pols,tmp,offset);
      for (size_t i(0);i<ress.size();++i)
        ress[i]+=std::complex<double>
          ((double)tmp[i].real(),(double)tmp[i].imag());
    }

  };// end of class Downcast_Contractor

  template <typename Derived>
  class Downcast_Contractor<Derived,double> {};

  template <class SType> inline void 
  Current::Contract(const Current &c,const Int_Vector &pols,
		    std::vector<std::complex<SType> > &ress,
		    const size_t &offset) const
  { dynamic_cast<const Current_Contractor<SType>*>
      (this)->SContract(c,pols,ress,offset); }

  typedef ATOOLS::Getter_Function<Current,Current_Key,
				  std::less<std::string> > Current_Getter;

  std::ostream &operator<<(std::ostream &str,const Current &c);

  typedef std::vector<Current*> Current_Vector;

}// end of namespace METOOLS

#endif
