#ifndef METOOLS_Explicit_Lorentz_Calculator_H
#define METOOLS_Explicit_Lorentz_Calculator_H

#include "METOOLS/Explicit/Color_Calculator.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/Vector.H"

namespace METOOLS {

  class Lorentz_Calculator {
  protected:

    Vertex *p_v;

    std::vector<int> m_r;

  public:

    // constructor
    Lorentz_Calculator(const Vertex_Key &key);

    // destructor
    virtual ~Lorentz_Calculator();

    // member functions
    virtual std::string Label() const = 0;

    virtual CObject *Evaluate(const CObject_Vector &j) = 0;

    virtual void Evaluate();

  };// end of class Lorentz_Calculator

  typedef std::vector<Lorentz_Calculator*> LC_Vector;

  typedef ATOOLS::Getter_Function<Lorentz_Calculator,Vertex_Key,
				  std::less<std::string> > LC_Getter;

}// end of namespace METOOLS

#endif
