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
  protected:

    Vertex *p_v;

    Complex m_cpl;

    int m_stat;

    static size_t s_cimin, s_cimax;

  public:

    // constructor
    inline Color_Calculator(const Vertex_Key &key): 
      p_v(key.p_v), m_cpl(1.0), m_stat(1) {}

    // destructor
    virtual ~Color_Calculator();

    // member functions
    virtual std::string Label() const = 0;

    virtual bool Evaluate(const CObject_Vector &j);
    virtual void AddJ(CObject *const c) = 0;

    inline static void SetCIMin(const size_t &cimin) { s_cimin=cimin; }
    inline static void SetCIMax(const size_t &cimax) { s_cimax=cimax; }

    inline static size_t CIMin() { return s_cimin; }
    inline static size_t CIMax() { return s_cimax; }

    inline Complex Coupling() const { return m_cpl; }

    inline int Stat() const { return m_stat; }

  };// end of class Color_Calculator

  typedef std::vector<Color_Calculator*> CC_Vector;

  typedef ATOOLS::Getter_Function<Color_Calculator,Vertex_Key,
				  std::less<std::string> > CC_Getter;

}// end of namespace METOOLS

#endif
