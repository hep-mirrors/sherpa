#ifndef Axis_C
#define Axis_C

#include "Axis.H"

namespace ATOOLS {

  template <class Value_Type>
  Axis<Value_Type>::Axis():
    m_variable(ATOOLS::Variable()),
    m_scalingmode(Reference),
    p_scaling(new ATOOLS::Id_Scaling<ValueType>()) {}

  template <class Value_Type>
  Axis<Value_Type>::~Axis()
  { delete p_scaling; }
	      
  template <class Value_Type>
  Value_Type Axis<Value_Type>::DisplayedValue(ValueType realvalue,ScalingMode tempsmode)
  {
    if (tempsmode==Unknown) tempsmode=m_scalingmode;
    switch (tempsmode) {
    case Reference:
      return (*p_scaling)(realvalue);
      break;
    case Unknown:
    case Identical:
      return realvalue;
      break;
    }
    return (ValueType)0.0;
  }
  
  template <class Value_Type>
  Value_Type Axis<Value_Type>::RealValue(ValueType displayedvalue,ScalingMode tempsmode)
  {
    if (tempsmode==Unknown) tempsmode=m_scalingmode;
    switch (tempsmode) {
    case Reference:
      return (*p_scaling)[displayedvalue];
      break;
    case Unknown:
    case Identical:
      return displayedvalue;
      break;
    }
    return (ValueType)0.0;
  }

} // end of namespace ATOOLS

#endif
