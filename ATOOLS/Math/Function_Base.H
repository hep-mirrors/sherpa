#ifndef ATOOLS_Math_Function_Base_H
#define ATOOLS_Math_Function_Base_H

#include <string>

namespace ATOOLS {

  class Function;

  class Function_Base {
  protected :

    double      m_defval;

    std::string m_type, m_name;

  public:

    // destructor
    virtual ~Function_Base();

    // setters
    virtual void SetParameters(double *parameters);

    Function *GetAIFunction();
    Function *GetAIGMeanFunction();

    // getters
    virtual std::string Type();

    // member functions
    virtual double GetValue(double x);
    virtual double operator()(double x);
    virtual double operator()();

    double FindZero(double, double, int MAX_ITR=150, double precision=2.2e-16);
    double WDBSolve(const double &y,const double &xmin,const double &xmax,
		    const double &precision=1.0e-12,const int maxit=100);

    inline void SetType(const std::string &type) { m_type=type; }

    inline std::string Name() const { return m_name; }

    inline void SetDefault(const double &defval) { m_defval=defval; }

    inline double Default() const { return m_defval; }
    
  };

// --------------------------------------------------
//           Doxygen part starts here
// --------------------------------------------------


/*!
 \file
 \brief contains the class Function_Base
*/


/*!
 \class Function_Base
 \brief this is a pure virtual base class for any function \f$y=f(x)\f$.

 This class provides an common interface to all function that 
 are only dependent on one parameter. All classes derived
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
