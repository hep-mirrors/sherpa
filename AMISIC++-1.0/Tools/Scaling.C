#ifndef Scaling_C
#define Scaling_C

#include "Scaling.H"
#include "Message.H"

namespace ATOOLS {

  template <class Value_Type>
  Value_Type Scaling_Base<Value_Type>::operator()(ValueType x)
  { 
    ATOOLS::msg.Error()<<"Scaling_Base::operator(): "
		       <<"Virtual method called!"<<std::endl; 
    return (Value_Type)0.0;
  }
  
  template <class Value_Type>
  Value_Type Scaling_Base<Value_Type>::operator[](ValueType x)
  { 
    ATOOLS::msg.Error()<<"Scaling_Base::operator[]: "
		       <<"Virtual method called!"<<std::endl; 
    return (Value_Type)0.0;
  }
  
  template <class Value_Type>
  const std::string Scaling_Base<Value_Type>::Name() const
  { 
    ATOOLS::msg.Error()<<"Scaling_Base::Name(): "
		       <<"Virtual method called!"<<std::endl; 
    return std::string("Scaling_Base");
  }
  
  template <class Value_Type>
  const std::string Id_Scaling<Value_Type>::Name() const 
  { return std::string("Id"); }

  template <class Value_Type>
  Value_Type Id_Scaling<Value_Type>::operator()(Value_Type x) 
  { return x; }

  template <class Value_Type>
  Value_Type Id_Scaling<Value_Type>::operator[](Value_Type y) 
  { return y; }

  template <class Value_Type>
  const std::string Log_Scaling<Value_Type>::Name() const 
  { return std::string("Log"); }

  template <class Value_Type>
  Value_Type Log_Scaling<Value_Type>::operator()(Value_Type x) 
  { return log(x); }

  template <class Value_Type>
  Value_Type Log_Scaling<Value_Type>::operator[](Value_Type y) 
  { return exp(y); }
    
  template <class Value_Type>
  const std::string Exp_Scaling<Value_Type>::Name() const 
  { return std::string("Exp"); }

  template <class Value_Type>
  Value_Type Exp_Scaling<Value_Type>::operator()(Value_Type x) 
  { return exp(x); }
    
  template <class Value_Type>
  Value_Type Exp_Scaling<Value_Type>::operator[](Value_Type y) 
  { return log(y); }

  template <class Value_Type>
  const std::string Sqr_Scaling<Value_Type>::Name() const 
  { return std::string("Sqr"); }

  template <class Value_Type>
  Value_Type Sqr_Scaling<Value_Type>::operator()(Value_Type x) 
  { return x*x; }

  template <class Value_Type>
  Value_Type Sqr_Scaling<Value_Type>::operator[](Value_Type y) 
  { return sqrt(y); }
    
  template <class Value_Type>
  const std::string Sqrt_Scaling<Value_Type>::Name() const 
  { return std::string("Sqrt"); }

  template <class Value_Type>
  Value_Type Sqrt_Scaling<Value_Type>::operator()(Value_Type x) 
  { return sqrt(x); }

  template <class Value_Type>
  Value_Type Sqrt_Scaling<Value_Type>::operator[](Value_Type y) 
  { return y*y; }

  template <class Value_Type>
  Log_B_Scaling<Value_Type>::Log_B_Scaling(Value_Type _b) 
  { b=_b; logb=log(b); }

  template <class Value_Type>
  const std::string Log_B_Scaling<Value_Type>::Name() const 
  { return std::string("Log_B_")+ToString(b); }

  template <class Value_Type>
  Value_Type Log_B_Scaling<Value_Type>::operator()(Value_Type x) 
  { return log(x)/logb; }

  template <class Value_Type>
  Value_Type Log_B_Scaling<Value_Type>::operator[](Value_Type y) 
  { return pow(b,y); }
    
  template <class Value_Type>
  B_To_X_Scaling<Value_Type>::B_To_X_Scaling(Value_Type _b) 
  { b=_b; }

  template <class Value_Type>
  const std::string B_To_X_Scaling<Value_Type>::Name() const 
  { return std::string("B_To_X_")+ToString(b); }

  template <class Value_Type>
  Value_Type B_To_X_Scaling<Value_Type>::operator()(Value_Type x) 
  { return pow(b,x); }

  template <class Value_Type>
  Value_Type B_To_X_Scaling<Value_Type>::operator[](Value_Type y) 
  { return log(y)/log(b); }
    
  template <class Value_Type>
  X_To_P_Scaling<Value_Type>::X_To_P_Scaling(Value_Type _p) { p=_p; }

  template <class Value_Type>
  const std::string X_To_P_Scaling<Value_Type>::Name() const 
  { return std::string("X_To_P_")+ToString(p); }

  template <class Value_Type>
  Value_Type X_To_P_Scaling<Value_Type>::operator()(Value_Type x) 
  { return pow(x,p); }

  template <class Value_Type>
  Value_Type X_To_P_Scaling<Value_Type>::operator[](Value_Type y) 
  { return pow(y,(Value_Type)1.0/p); }

}

#endif
