#ifndef Function_Base_H
#define Function_Base_H

namespace AMATOOLS {
  class Function_Base {
  public:
    virtual void   SetParameters(double *) { };
    virtual double GetValue(double xx) {return (*this)(xx);}   
    virtual double operator() (double)=0;
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
