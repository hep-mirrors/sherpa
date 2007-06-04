#ifndef Vertex_H
#define Vertex_H

#include "Current_Base.H"

namespace EXTRAXS {

  class Vertex {
  protected:

    Current_Base *p_a, *p_b, *p_c;

    bool   m_zero, m_sign;
    size_t m_fperm;

    template <typename CType> inline void 
    AddJ(const CType &c) const { p_c->AddJ(m_sign?-c:c); }

    virtual std::string VLabel() const;

  public:

    Vertex(Current_Base *const c);

    virtual ~Vertex();

    virtual void Evaluate() = 0;

    void FindPermutation();

    std::string VId() const;

    void CollectGraphs(Graph_Node *graph) const;

    inline void SetJA(Current_Base *const a) { p_a=a; a->AttachOut(this); }
    inline void SetJB(Current_Base *const b) { p_b=b; b->AttachOut(this); }
    inline void SetJC(Current_Base *const c) { p_c=c; c->AttachIn(this);  }

    inline void SetZero(const bool zero) { m_zero=zero; }

    inline Current_Base *const JA() const { return p_a; }
    inline Current_Base *const JB() const { return p_b; }
    inline Current_Base *const JC() const { return p_c; }

    inline bool Zero() const { return m_zero; }
    inline bool Sign() const { return m_sign; }

    inline size_t FPerm() const { return m_fperm; }

  };// end of class Vertex

  struct Vertex_Key {
    Current_Base  *p_a, *p_b, *p_c;
    std::string    m_model;
    XS_Model_Base *p_model;
    inline Vertex_Key(Current_Base *const a,Current_Base *const b,
		      Current_Base *const c,const std::string &cmodel,
		      XS_Model_Base *const model): 
      p_a(a), p_b(b), p_c(c), m_model(cmodel), p_model(model) {}
    std::string Type() const;
    inline Vertex_Key SwapAB() 
    { return Vertex_Key(p_b,p_a,p_c,m_model,p_model); }
    bool operator<(const Vertex_Key &k) const;
  };// end of struct Vertex_Key

  typedef ATOOLS::Getter_Function<Vertex,Vertex_Key> Vertex_Getter;

  std::ostream &operator<<(std::ostream &str,const Vertex &v);

}// end of namespace EXTRAXS

#endif
