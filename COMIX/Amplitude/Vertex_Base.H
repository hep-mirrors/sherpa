#ifndef COMIX_Amplitude_Vertex_Base_H
#define COMIX_Amplitude_Vertex_Base_H

#include "COMIX/Amplitude/Current_Base.H"

namespace COMIX {

  struct Vertex_Key {
    Current_Base  *p_a, *p_b, *p_c;
    Model *p_model;
    inline Vertex_Key(Current_Base *const a,Current_Base *const b,
		      Current_Base *const c,Model *const model): 
      p_a(a), p_b(b), p_c(c), p_model(model) {}
    std::string Type() const;
    std::string ID(const int dir=0) const;
    inline Vertex_Key SwapAB() 
    { return Vertex_Key(p_b,p_a,p_c,p_model); }
    bool operator<(const Vertex_Key &k) const;
  };// end of struct Vertex_Key

  class Vertex_Base {
  protected:

    Current_Base *p_a, *p_b, *p_c;

    bool   m_sign, m_act, m_zero;
    size_t m_fperm, m_oew, m_oqcd;

    std::string m_tag;

    double m_cplfac;

    static size_t s_vlmode, s_cimin, s_cimax;

    template <typename CType> inline void 
    AddJ(const CType &c) const { p_c->AddJ(m_sign?-c:c); }

    virtual std::string CVLabel() const;

  public:

    // constructor
    Vertex_Base(const Vertex_Key &key,
		const size_t &oew=0,const size_t &oqcd=0);

    // destructor
    virtual ~Vertex_Base();

    // member functions
    virtual void Evaluate() = 0;

    void FindPermutation();

    bool Map(const Vertex_Base &v);

    std::string VId() const;
    std::string VLabel() const;

    void CollectGraphs(Graph_Node *graph) const;

    // inline functions
    inline void SetCplFac(const double &fac) { m_cplfac=fac; }

    inline void SetJA(Current_Base *const a) 
    { p_a=a; if (a!=NULL) a->AttachOut(this); }
    inline void SetJB(Current_Base *const b) 
    { p_b=b; if (b!=NULL) b->AttachOut(this); }
    inline void SetJC(Current_Base *const c) 
    { p_c=c; if (c!=NULL) c->AttachIn(this); }

    inline Current_Base *const JA() const { return p_a; }
    inline Current_Base *const JB() const { return p_b; }
    inline Current_Base *const JC() const { return p_c; }

    inline bool Zero() const { return m_zero; }
    inline void SetZero()    { m_zero=true;   }

    inline bool Sign() const { return m_sign; }

    inline bool Active() const { return m_act; }

    inline size_t FPerm() const { return m_fperm; }

    inline size_t OrderEW() const  { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }

    inline std::string Tag() const { return m_tag; }

    inline static void SetVLMode(const size_t &mode) { s_vlmode=mode; }

    inline static size_t VLMode() { return s_vlmode; }

    inline static void SetCIMin(const size_t &cimin) { s_cimin=cimin; }
    inline static void SetCIMax(const size_t &cimax) { s_cimax=cimax; }

  };// end of class Vertex_Base
  
  typedef ATOOLS::Getter_Function<Vertex_Base,Vertex_Key,
				  std::less<std::string> > Vertex_Getter;
  typedef ATOOLS::Getter_Function<void,const Model*,
				  std::less<std::string> > Filler_Getter;

  std::ostream &operator<<(std::ostream &str,const Vertex_Base &v);

}// end of namespace COMIX

#endif
