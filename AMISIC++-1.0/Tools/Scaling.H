#ifndef Scaling_H
#define Scaling_H

#include "My_String.H"
#include "My_Stringstream.H"
#include "My_Functional.H"

namespace ATOOLS {

  template <class Value_Type>
  std::string ToString(const Value_Type value) {
    std::stringstream converter;
    std::string converted;
    converter<<value;
    converter>>converted;
    return converted;
  }

  template <class Value_Type>
  class Scaling_Base: public std::unary_function<Value_Type,Value_Type> {
  public:

    typedef Value_Type ValueType;

  public:

    // getter
    virtual const std::string Name() const;

    // member functions
    virtual Value_Type operator()(Value_Type x);
    virtual Value_Type operator[](Value_Type y);

  };
  
  template <class Value_Type>
  class Id_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const { return std::string("Id"); }

    Value_Type operator()(Value_Type x) { return x; }
    Value_Type operator[](Value_Type y) { return y; }
    
  };  

  template <class Value_Type>
  class Log_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const { return std::string("Log"); }

    Value_Type operator()(Value_Type x) { return log(x); }
    Value_Type operator[](Value_Type y) { return exp(y); }
    
  };  

  template <class Value_Type>
  class Sqr_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const { return std::string("Sqr"); }

    Value_Type operator()(Value_Type x) { return x*x; }
    Value_Type operator[](Value_Type y) { return sqrt(y); }
    
  };  

  template <class Value_Type>
  class Log_B_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b;

  public:

    Log_B_Scaling(Value_Type _b) { b=_b; }

    const std::string Name() const { return std::string("Log_B_")+ToString(b); }

    Value_Type operator()(Value_Type x) { return log(x)/log(b); }
    Value_Type operator[](Value_Type y) { return pow(b,y); }
    
  };  

  template <class Value_Type>
  class B_To_X_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b;

  public:

    B_To_X_Scaling(Value_Type _b) { b=_b; }

    const std::string Name() const { return std::string("B_To_X_")+ToString(b); }

    Value_Type operator()(Value_Type x) { return pow(b,x); }
    Value_Type operator[](Value_Type y) { return log(y)/log(b); }
    
  };  

  template <class Value_Type>
  class X_To_P_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type p;

  public:

    X_To_P_Scaling(Value_Type _p) { p=_p; }

    const std::string Name() const { return std::string("X_To_P_")+ToString(p); }

    Value_Type operator()(Value_Type x) { return pow(x,p); }
    Value_Type operator[](Value_Type y) { return pow(y,(Value_Type)1.0/p); }
    
  };  

  //////////////////
  // Doxygen part //
  //////////////////

  /*
    \class Scaling_Base
    \brief Respresents a scaling function

    The classes derived from Scaling_Base are employed for the scaling of values 
    within an Axis object. There exist different implementations of scaling functions.
    Namely there are

    \itemize
      \item identity  (x'=x)
      \item square    (x'=x^2)
      \item power     (x'=x^p)
      \item logarithm (x'=log(x))
      \item logarithm (x'=log_{b}(x))
    \enditemize
    
    and their inverse functions respectively.
    The classes derived from Scaling_Base own the inverse scaling functions as well.
  */

  /*
    \fn Scaling_Base::Name()
    \brief Returns the name of the Scaling_Base
  */

  /*
    \fn Scaling_Base::operator()(Value_Type x)
    \brief Performs the scaling of x.
  */

  /*
    \fn Scaling_Base::operator()(Value_Type x)
    \brief Performs the inverse scaling of x.
  */

  /*
    \class Id_Scaling
    \brief Provides identical scaling.
  */

  /*
    \class Log_Scaling
    \brief Provides logarithmic scaling.
  */

  /*
    \class Sqr_Scaling
    \brief Provides scaling by x^2.
  */

  /*
    \class Log_B_Scaling
    \brief Provides scaling by \log_b(x).
  */

  /*
    \class B_To_X_Scaling
    \brief Provides scaling by b^x.
  */

  /*
    \class X_To_P_Scaling
    \brief Provides scaling by x^p.
  */

} // end of namespace ATOOLS

#include "Scaling.C"

#endif
