#ifndef Debugger_inl_H
#define Debugger_inl_H

#include "Debugger.H"
#include "Type.H"

namespace ATOOLS {

  template <class ValueType>
  Value_Carrier<ValueType> *Value_Carrier_Base::UpCast()
  {
    return dynamic_cast<Value_Carrier<ValueType>*>(this);
  }

  template <class ValueType> 
  ValueType *Value_Carrier_Base::Address() 
  { 
    return UpCast<ValueType>()->Address(); 
  }

  template <class ValueType> 
  ValueType Value_Carrier_Base::Value() 
  { 
    return UpCast<ValueType>()->Value(); 
  }

  template <class ValueType>
  Value_Carrier<ValueType>::Value_Carrier(Value_Type &address): 
    p_address(&address) 
  { 
    Value_Type value;
    m_type=Type::GetType(value);
  }
  
  template <class ValueType>
  ValueType *Value_Carrier<ValueType>::Address() const
  {
    return p_address;
  }
  
  template <class ValueType>
  ValueType Value_Carrier<ValueType>::Value() const
  {
    return *p_address;
  }  
  
  template <class ValueType> 
  void Value_Carrier<ValueType>::PrintAddress(std::ostream &ostr) const
  {
    ostr<<Address();
  }
  
  template <class ValueType> 
  void Value_Carrier<ValueType>::PrintValue(std::ostream &ostr) const
  {
    ostr<<Value();
  }
  
}// end of namespace ATOOLS

#endif
