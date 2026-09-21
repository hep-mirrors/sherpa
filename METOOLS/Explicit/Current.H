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
    // Double-double shadow of m_p, used only for the p^2 bookkeeping below.
    // m_p itself stays double on purpose.
    ATOOLS::Vec4<ATOOLS::DDouble> m_ph;
    // Externally supplied wide momentum (Amplitude's joint on-shell /
    // conserving projection). When set it wins over the per-leg
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

    inline void SetPWideExt(const ATOOLS::Vec4<ATOOLS::DDouble> &ph)
    { m_phx=ph; m_hasphx=true; }

    inline void SetPWide()
    {
      if (m_hasphx) { m_ph=m_phx; return; }
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

    inline double P2() const { return (double)m_p2; }

    inline const ATOOLS::Vec4<ATOOLS::DDouble> &PH() const { return m_ph; }

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
