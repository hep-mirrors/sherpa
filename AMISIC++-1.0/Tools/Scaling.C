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
  
}

#endif
