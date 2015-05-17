#ifndef METOOLS_Explicit_Vertex_H
#define METOOLS_Explicit_Vertex_H

#include "METOOLS/Explicit/Vertex_Key.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "METOOLS/Explicit/Lorentz_Calculator.H"

namespace METOOLS {

  typedef std::vector<Complex> Complex_Vector;
  typedef std::vector<Complex_Vector> Complex_Matrix;

  class Dipole_Info;
  class Dipole_Kinematics;

  class Vertex {
  protected:

    MODEL::Single_Vertex *p_v;

    Current_Vector  m_j;
    Current        *p_c;

    Dipole_Info       *p_info;
    Dipole_Kinematics *p_kin;

    LC_Vector m_lc;
    CC_Vector m_cc;

    SizeT_Vector m_h, m_hjc, m_cjc;

    bool   m_sign, m_act, m_zero;
    size_t m_fperm;

    std::vector<int> m_order;

    double m_icplfac;

    std::vector<const CObject_Vector*> m_hjj;
    CObject_Vector m_cjj;

    static size_t s_vlmode;

    std::string CVLabel() const;

    friend class Lorentz;
    friend class Color;

  public:

    // constructor
    Vertex(const Vertex_Key &sv);

    // destructor
    ~Vertex();

    // member functions
    void Evaluate();

    void FindPermutation();
    void InitPols();

    bool Map(const Vertex &v);

    std::string VId() const;
    std::string VLabel() const;

    void CollectGraphs(Graph_Node *graph) const;

    void AddJ(const Current_Vector &j);

    // inline functions
    inline void AddJ(Current *const j) 
    { if (j) { m_j.push_back(j); j->AttachOut(this); } }
    inline Current *J(const size_t &i) const { return m_j[i]; }
    inline const Current_Vector &J() const { return m_j; }

    inline void SetJC(Current *const c,const bool a=true) 
    { p_c=c; if (c!=NULL && a) c->AttachIn(this); }
    inline Current *JC() const { return p_c; }

    inline void ClearJ() { m_j.clear(); }

    inline Dipole_Info       *Info() const { return p_info; }
    inline Dipole_Kinematics *Kin() const  { return p_kin;  }

    inline void AddJ(CObject *const c) const 
    {
      if (m_sign) c->Invert();
      if (m_icplfac!=1.0) c->Divide(m_icplfac);
      p_c->AddJ(c);
    }

    inline void SetCplFac(const double &fac) { m_icplfac=1.0/fac; }

    inline void SetZero(const bool zero=true) { m_zero=zero; }

    inline bool Zero() const { return m_zero; }

    inline bool Sign() const { return m_sign; }

    inline bool Active() const { return m_act; }

    inline size_t FPerm() const { return m_fperm; }

    inline const std::vector<int> &Order() const { return m_order; }
    inline int Order(const size_t &id) const { return m_order[id]; }

    inline static void SetVLMode(const size_t &mode) { s_vlmode=mode; }

    inline static size_t VLMode() { return s_vlmode; }

    inline const LC_Vector &Lorentz() const { return m_lc; }
    inline const CC_Vector &Color() const   { return m_cc; }

    inline size_t H(const size_t &i) const { return m_h[i]; }

    inline MODEL::Single_Vertex *V() const { return p_v; }

  };// end of class Vertex
  
  std::ostream &operator<<(std::ostream &str,const Vertex &v);

}// end of namespace METOOLS

#endif
