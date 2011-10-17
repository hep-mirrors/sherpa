#ifndef METOOLS_Explicit_Vertex_H
#define METOOLS_Explicit_Vertex_H

#include "METOOLS/Explicit/Vertex_Key.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "METOOLS/Explicit/Lorentz_Calculator.H"

namespace METOOLS {

  typedef std::vector<Complex> Complex_Vector;

  class Dipole_Info;
  class Dipole_Kinematics;

  class Vertex {
  protected:

    Current *p_a, *p_b, *p_c, *p_e;

    Dipole_Info       *p_info;
    Dipole_Kinematics *p_kin;

    LC_Vector m_lc;
    CC_Vector m_cc;

    Complex_Vector m_cpl;
    SizeT_Vector   m_h;

    bool   m_sign, m_act, m_zero;
    size_t m_fperm, m_oew, m_oqcd;

    double m_icplfac;

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

    Current *J(const size_t &i) const;

    // inline functions
    inline void SetJA(Current *const a) 
    { p_a=a; if (a!=NULL) a->AttachOut(this); }
    inline void SetJB(Current *const b) 
    { p_b=b; if (b!=NULL) b->AttachOut(this); }
    inline void SetJE(Current *const e) 
    { p_e=e; if (e!=NULL) e->AttachOut(this); }
    inline void SetJC(Current *const c,const bool a=true) 
    { p_c=c; if (c!=NULL && a) c->AttachIn(this); }

    inline Current *JA() const { return p_a; }
    inline Current *JB() const { return p_b; }
    inline Current *JE() const { return p_e; }
    inline Current *JC() const { return p_c; }

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

    inline size_t OrderEW() const  { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }

    inline static void SetVLMode(const size_t &mode) { s_vlmode=mode; }

    inline static size_t VLMode() { return s_vlmode; }

    inline const LC_Vector &Lorentz() const { return m_lc; }
    inline const CC_Vector &Color() const   { return m_cc; }

    inline const Complex &Coupling(const size_t &i) const 
    { return m_cpl[i]; }

    inline size_t H(const size_t &i) const { return m_h[i]; }

  };// end of class Vertex
  
  std::ostream &operator<<(std::ostream &str,const Vertex &v);

}// end of namespace METOOLS

#endif
