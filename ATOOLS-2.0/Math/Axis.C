#include "Axis.H"

#include "Data_Reader.H"

using namespace ATOOLS;

template <class ValueType>
Axis<ValueType>::Axis():
  m_scalingmode(Reference),
  p_variable(Variable_Getter::GetObject("","")),
  p_scaling(Scaling_Getter::GetObject("","")) {}

template <class ValueType>
Axis<ValueType>::Axis(const Axis &ref):
  m_scalingmode(ref.m_scalingmode),
  p_variable(Variable_Getter::GetObject(ref.p_variable->Name(),
					ref.p_variable->Name())),
  p_scaling(Scaling_Getter::GetObject(ref.p_scaling->Name(),
				      ref.p_scaling->Name())) {}

template <class ValueType>
Axis<ValueType>::~Axis() {}

template <class ValueType>
ValueType Axis<ValueType>::DisplayedValue(const Value_Type &realvalue,
					  ScalingModeID tempsmode)
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
  return (Value_Type)0.0;
}

template <class ValueType>
ValueType Axis<ValueType>::RealValue(const Value_Type &displayedvalue,
				     ScalingModeID tempsmode)
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
  return (Value_Type)0.0;
}

template <class ValueType>
void Axis<ValueType>::SetScaling(const std::string &scalename)
{
  p_scaling=Scaling_Getter::GetObject(scalename,scalename);
  if (p_scaling==NULL) p_scaling=Scaling_Getter::GetObject("","");
}

template <class ValueType>
void Axis<ValueType>::SetVariable(const std::string &variablename)
{
  p_variable=Variable_Getter::GetObject(variablename,variablename);
  if (p_variable==NULL) p_variable=Variable_Getter::GetObject("","");
}

template class Axis<double>;
