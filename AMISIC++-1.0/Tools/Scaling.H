#ifndef Scaling_H
#define Scaling_H

#include "My_String.H"
#include "My_Stringstream.H"
#include "My_Functional.H"

namespace ATOOLS {

  template <class Value_Type>
  std::string ToString(Value_Type value) {
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
    virtual std::string Name();

    // member functions
    virtual Value_Type operator()(Value_Type x);
    virtual Value_Type operator[](Value_Type y);

  };
  
  template <class Value_Type>
  class Id_Scaling: public Scaling_Base<Value_Type> {
  public:

    std::string Name() { return std::string("Id"); }

    Value_Type operator()(Value_Type x) { return x; }
    Value_Type operator[](Value_Type y) { return y; }
    
  };  

  template <class Value_Type>
  class Log_Scaling: public Scaling_Base<Value_Type> {
  public:

    std::string Name() { return std::string("Log"); }

    Value_Type operator()(Value_Type x) { return log(x); }
    Value_Type operator[](Value_Type y) { return exp(y); }
    
  };  

  template <class Value_Type>
  class Sqr_Scaling: public Scaling_Base<Value_Type> {
  public:

    std::string Name() { return std::string("Sqr"); }

    Value_Type operator()(Value_Type x) { return x*x; }
    Value_Type operator[](Value_Type y) { return sqrt(y); }
    
  };  

  template <class Value_Type>
  class Log_B_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b;

  public:

    Log_B_Scaling(Value_Type _b) { b=_b; }

    std::string Name() { return std::string("Log_B_")+ToString(b); }

    Value_Type operator()(Value_Type x) { return log(x)/log(b); }
    Value_Type operator[](Value_Type y) { return pow(b,y); }
    
  };  

  template <class Value_Type>
  class B_To_X_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b;

  public:

    B_To_X_Scaling(Value_Type _b) { b=_b; }

    std::string Name() { return std::string("B_To_X_")+ToString(b); }

    Value_Type operator()(Value_Type x) { return pow(b,x); }
    Value_Type operator[](Value_Type y) { return log(y)/log(b); }
    
  };  

  template <class Value_Type>
  class X_To_P_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type p;

  public:

    X_To_P_Scaling(Value_Type _p) { p=_p; }

    std::string Name() { return std::string("X_To_P_")+ToString(p); }

    Value_Type operator()(Value_Type x) { return pow(x,p); }
    Value_Type operator[](Value_Type y) { return pow(y,(Value_Type)1.0/p); }
    
  };  

} // end of namespace ATOOLS

#include "Scaling.C"

#endif
