#ifndef Debugger_inl_H
#define Debugger_inl_H

#include "Debugger.H"

namespace ATOOLS {

  template <typename ValueType>
  Value_Carrier<ValueType> *Value_Carrier_Base::UpCast()
  {
    return dynamic_cast<Value_Carrier<ValueType>*>(this);
  }

  template <typename ValueType> 
  ValueType *Value_Carrier_Base::Address() 
  { 
    return UpCast<ValueType>()->Address(); 
  }

  template <typename ValueType> 
  ValueType Value_Carrier_Base::Value() 
  { 
    return UpCast<ValueType>()->Value(); 
  }

  template <typename ValueType>
  Value_Carrier<ValueType>::Value_Carrier(Value_Type &address): 
    p_address(&address) 
  { 
    m_type=typeid(address).name();
  }
  
  template <typename ValueType>
  ValueType *Value_Carrier<ValueType>::Address() const
  {
    return p_address;
  }
  
  template <typename ValueType>
  ValueType Value_Carrier<ValueType>::Value() const
  {
    return *p_address;
  }  
  
  template <typename ValueType> 
  void Value_Carrier<ValueType>::PrintAddress(std::ostream &ostr) const
  {
    ostr<<Address();
  }
  
  template <typename ValueType> 
  void Value_Carrier<ValueType>::PrintValue(std::ostream &ostr) const
  {
    ostr<<Value();
  }
  
  template <typename ValueType>
  void Debugger::AddValue(const std::string &spec,ValueType &value) 
  { 
    m_values[spec] = new Value_Carrier<ValueType>(value); 
  }

}// end of namespace ATOOLS

#endif
