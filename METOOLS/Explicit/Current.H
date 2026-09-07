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
    // Virtuality of m_p, accumulated recursively rather than recomputed as
    // m_p.Abs2(). For an internal line built from nearly-collinear externals
    // the two differ enormously: m_p.Abs2() forms an O(E^2) difference to
    // reach an O(E^2 theta^2) answer, so for the ISR photons of a lepton
    // collider (theta ~ m_e/E_beam) every significant digit of the propagator
    // denominator is lost. Evaluate() instead builds
    //   p^2 = sum_i p_i^2 + 2 sum_{i<j} p_i.p_j
    // from the sub-current virtualities, which are exact on-shell masses at
    // the leaves, so the only cancellation left sits inside a single
    // Vec4D::LCDot. See AddPropagator() in Currents/{F,V,S}_C.C.
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
    inline void SetP(const ATOOLS::Vec4D &p)
    { m_p=p; m_ph=p; m_p2=ATOOLS::sqr(m_mass); }

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
