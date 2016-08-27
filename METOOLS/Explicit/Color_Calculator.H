#ifndef METOOLS_Explicit_Color_Calculator_H
#define METOOLS_Explicit_Color_Calculator_H

#include "METOOLS/Explicit/Vertex_Key.H"
#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace METOOLS {

  class Vertex;
  class Current;

  class Color_Calculator {
  public:

    struct CInfo {
      int m_cr, m_ca;
      Complex m_s;
      inline CInfo(const int cr,const int ca,
		   const Complex &s=Complex(1.0)):
	m_cr(cr), m_ca(ca), m_s(s) {}
      inline int &operator()(const int i)
      { return i?m_ca:m_cr; }
      inline const int &operator()(const int i) const
      { return i?m_ca:m_cr; }
    };// end of struct CInfo

  protected:

    Vertex *p_v;
    Complex m_cpl;

    std::vector<CInfo> m_c;

    static size_t s_cimin, s_cimax;

  public:

    // constructor
    inline Color_Calculator(const Vertex_Key &key): 
      p_v(key.p_v), m_cpl(1.0) {}

    // destructor
    virtual ~Color_Calculator();

    virtual void AddJ(CObject *const j);

    // member functions
    virtual std::string Label() const = 0;

    virtual bool Evaluate(const CObject_Vector &j);

    inline static void SetCIMin(const size_t &cimin) { s_cimin=cimin; }
    inline static void SetCIMax(const size_t &cimax) { s_cimax=cimax; }

    inline static size_t CIMin() { return s_cimin; }
    inline static size_t CIMax() { return s_cimax; }

    inline Complex Coupling() const { return m_cpl; }

    inline int Stat() const { return m_c.size(); }

  };// end of class Color_Calculator

  std::ostream &operator<<(std::ostream &s,const Color_Calculator::CInfo &i);

  typedef std::vector<Color_Calculator*> CC_Vector;

  typedef ATOOLS::Getter_Function<Color_Calculator,Vertex_Key,
				  std::less<std::string> > CC_Getter;

}// end of namespace METOOLS

#endif
