#include "Variable.H"

#include "Message.H"
#include "Vector.H"

using namespace ATOOLS;
  
template <class ValueType>
Variable_Base<ValueType>::Variable_Base(const std::string &name):
  m_name(name),
  m_selectorid(-1) {}

template <class ValueType>
Variable_Base<ValueType>::~Variable_Base() {}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value(const Vec3D *vectors)
{
  msg.Error()<<"Variable_Base::Value("<<vectors<<"): "
	     <<"Virtual method called."<<std::endl;
  return 0.0;
}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value(const Vec4D *vectors)
{
  msg.Error()<<"Variable_Base::Value("<<vectors<<"): "
	     <<"Virtual method called."<<std::endl;
  return 0.0;
}

template <class ValueType>
void Variable_Base<ValueType>::ShowVariables(const int mode)
{
  if (!msg.LevelIsInfo() || mode==0) return;
  msg.Out()<<"Variable_Base::ShowVariables(): {\n\n";
  Variable_Getter::PrintGetterInfo(msg.Out(),20);
  msg.Out()<<"\n}"<<std::endl;
}

template <class ValueType>
class No_Variable: public Variable_Base<ValueType> {
public:
  No_Variable(): Variable_Base<ValueType>("") {}
};// end of class No_Variable
  
template <class ValueType>
class PPerp: public Variable_Base<ValueType> {
public:
  PPerp(): Variable_Base<ValueType>("p_\\perp") { m_selectorid=12; }
  ValueType Value(const Vec3D *vectors) 
  { return Vec4D(0.0,vectors[0]).PPerp(); }
  ValueType Value(const Vec4D *vectors) { return vectors[0].PPerp(); }
};// end of class PPerp
  
template <class ValueType>
class EPerp: public Variable_Base<ValueType> {
public:
  EPerp(): Variable_Base<ValueType>("E_\\perp") { m_selectorid=15; }
  ValueType Value(const Vec4D *vectors) { return vectors[0].EPerp(); }
};// end of class EPerp
  
template <class ValueType>
class Energy: public Variable_Base<ValueType> {
public:
  Energy(): Variable_Base<ValueType>("E") { m_selectorid=11; }
  ValueType Value(const Vec4D *vectors) { return vectors[0][0]; }
};// end of class Energy
  
template <class ValueType>
class Mass: public Variable_Base<ValueType> {
public:
  Mass(): Variable_Base<ValueType>("m") { m_selectorid=21; }
  ValueType Value(const Vec4D *vectors) { return vectors[0].Mass(); }
};// end of class Mass
  
template <class ValueType>
class Rapidity: public Variable_Base<ValueType> {
public:
  Rapidity(): Variable_Base<ValueType>("y") { m_selectorid=13; }
  ValueType Value(const Vec4D *vectors) { return vectors[0].Y(); }
};// end of class Rapidity
  
template <class ValueType>
class Eta: public Variable_Base<ValueType> {
public:
  Eta(): Variable_Base<ValueType>("\\eta") { m_selectorid=16; }
  ValueType Value(const Vec4D *vectors) { return vectors[0].Eta(); }
};// end of class Eta
  
template <class ValueType>
class Theta: public Variable_Base<ValueType> {
public:
  Theta(): Variable_Base<ValueType>("\\theta") { m_selectorid=14; }
  ValueType Value(const Vec3D *vectors) 
  { return Vec4D(0.0,vectors[0]).Theta(); }
  ValueType Value(const Vec4D *vectors) { return vectors[0].Theta(); }
};// end of class Theta
  
template <class ValueType>
class Theta2: public Variable_Base<ValueType> {
public:
  Theta2(): Variable_Base<ValueType>("\\theta_{ij}") { m_selectorid=22; }
  ValueType Value(const Vec3D *vectors) 
  { return Vec4D(0.0,vectors[1]).Theta(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors) 
  { return vectors[1].Theta(vectors[0]); }
};// end of class Theta2

template class Variable_Base<double>;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Variable_Base<double>
#define PARAMETER_TYPE std::string
#include "Getter_Function.C"

template <class Class>
Variable_Base<double> *GetVariable(const std::string &parameter) 
{
  return new Class();
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Variable_Base<double> *						\
  NAME::operator()(const std::string &parameter) const			\
  { return GetVariable<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<PRINT; }

#define DEFINE_VARIABLE_GETTER(CLASS,NAME,TAG,PRINT)			\
  template class CLASS;							\
  DECLARE_GETTER(NAME,TAG,Variable_Base<double>,std::string);		\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

DEFINE_VARIABLE_GETTER(No_Variable<double>,No_Variable_Getter,
		       "","")
DEFINE_VARIABLE_GETTER(PPerp<double>,PPerp_Getter,
		       "p_\\perp","p_\\perp")
DEFINE_VARIABLE_GETTER(EPerp<double>,EPerp_Getter,
		       "E_\\perp","E_\\perp")
DEFINE_VARIABLE_GETTER(Energy<double>,Energy_Getter,
		       "E","E")
DEFINE_VARIABLE_GETTER(Mass<double>,Mass_Getter,
		       "m","m")
DEFINE_VARIABLE_GETTER(Rapidity<double>,Rapidity_Getter,
		       "y","y")
DEFINE_VARIABLE_GETTER(Eta<double>,Eta_Getter,
		       "\\eta","\\eta")
DEFINE_VARIABLE_GETTER(Theta<double>,Theta_Getter,
		       "\\theta","\\theta")
DEFINE_VARIABLE_GETTER(Theta2<double>,Theta2_Getter,
		       "\\theta_{ij}","\\theta_{ij}")

