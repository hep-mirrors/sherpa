#ifndef Axis_C
#define Axis_C

#include "Axis.H"
#include "Data_Reader.H"

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
  Value_Type Axis<Value_Type>::DisplayedValue(ValueType realvalue,ScalingModeID tempsmode)
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
  Value_Type Axis<Value_Type>::RealValue(ValueType displayedvalue,ScalingModeID tempsmode)
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

  template <class Value_Type>
  void Axis<Value_Type>::SetScaling(std::string scalename)
  {
    ValueType argx;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
    reader->SetString(scalename);
    if (scalename==std::string("Id")) {
      SetScaling(new ATOOLS::Id_Scaling<ValueType>()); 
    }
    if (scalename==std::string("Log")) {
      SetScaling(new ATOOLS::Log_Scaling<ValueType>()); 
    }
    else if (scalename==std::string("Exp")) {
      SetScaling(new ATOOLS::Exp_Scaling<ValueType>()); 
    }
    else if (scalename==std::string("Sqr")) {
      SetScaling(new ATOOLS::Sqr_Scaling<ValueType>()); 
    }
    else if (scalename==std::string("Sqrt")) {
      SetScaling(new ATOOLS::Sqrt_Scaling<ValueType>()); 
    }
    else if (reader->ReadFromString(argx,"Log_B_")) { 
      SetScaling(new ATOOLS::Log_B_Scaling<ValueType>(argx)); 
    }
    else if (reader->ReadFromString(argx,"B_To_X_")) { 
      SetScaling(new ATOOLS::B_To_X_Scaling<ValueType>(argx)); 
    }
    else if (reader->ReadFromString(argx,"X_To_P_")) { 
      SetScaling(new ATOOLS::X_To_P_Scaling<ValueType>(argx)); 
    }
    delete reader;
  }

} // end of namespace ATOOLS

#endif
