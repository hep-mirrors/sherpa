#ifndef Function_Base_H
#define Function_Base_H

#include<string>

namespace AMATOOLS {
  class Function_Base {
  protected :
    double      m_defval;
    std::string m_type;
  public:
    virtual void        SetDefault(double _def)    { m_defval = _def; }
    virtual void        SetType(std::string _type) { m_type   = _type; }
    virtual void        SetParameters(double *)    { }
    virtual double      GetValue(double xx)        { return (*this)(xx); }   
    virtual double      operator()(double)         { return m_defval; }
    virtual double      operator()()               { return m_defval; }
    virtual std::string Type()                     { return m_type; }
  };

// --------------------------------------------------
//           Doxygen part starts here
// --------------------------------------------------


/*!
 \file
 \brief contains the class AMATOOLS::Function_Base
*/


/*!
 \class Function_Base
 \brief this is a pure virtual base class for any function \f$y=f(x)\f$.

 This class provides an common interface to all function that 
 are only dipendent on one parameter. All classes derived
 from this class might e.g. be itegrated by Gauss_Integrater.
*/

/*!
 \fn virtual void Function_Base::SetParameters(double *)
  \brief  possibility to pass a number of paremters to the function.
*/

/*!
 \fn virtual double Function_Base::GetValue(double xx)
 \brief synonym for calling Function_Base::operator()();
*/

/*!
 \fn virtual double Function_Base::operator() (double x)
 \brief returns a value for \f$f(x)\f$
*/

}

#endif
