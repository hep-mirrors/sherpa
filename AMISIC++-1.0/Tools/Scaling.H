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

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class Log_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class Exp_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class Sqr_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class Sqrt_Scaling: public Scaling_Base<Value_Type> {
  public:

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class Log_B_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b, logb;

  public:

    Log_B_Scaling(Value_Type _b);

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class B_To_X_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type b;

  public:

    B_To_X_Scaling(Value_Type _b);

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

  template <class Value_Type>
  class X_To_P_Scaling: public Scaling_Base<Value_Type> {
  private:
    
    Value_Type p;

  public:

    X_To_P_Scaling(Value_Type _p);

    const std::string Name() const;

    Value_Type operator()(Value_Type x);
    Value_Type operator[](Value_Type y);
    
  };  

} // end of namespace ATOOLS

#include "Scaling.C"

#endif
