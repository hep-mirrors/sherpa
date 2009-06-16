#include "ATOOLS/Math/Scaling.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS; 

template <class ValueType>
Scaling_Base<ValueType>::~Scaling_Base() {}

template <class ValueType>
ValueType Scaling_Base<ValueType>::operator()(const Value_Type &x)
{ 
  msg_Error()<<"Scaling_Base::operator(): "
		       <<"Virtual method called!"<<std::endl; 
  return (Value_Type)0.0;
}

template <class ValueType>
ValueType Scaling_Base<ValueType>::operator[](const Value_Type &y)
{ 
  msg_Error()<<"Scaling_Base::operator[]: "
		     <<"Virtual method called!"<<std::endl; 
  return (Value_Type)0.0;
}

template <class ValueType>
const std::string Scaling_Base<ValueType>::Name() const 
{
  return m_name; 
}

template <class ValueType>
void Scaling_Base<ValueType>::ShowScalings(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<"Scaling_Base::ShowScalings(): {\n\n";
  Scaling_Getter::PrintGetterInfo(msg->Out(),20);
  msg_Out()<<"\n}"<<std::endl;
}

template <class Value_Type>
class Id_Scaling: public Scaling_Base<Value_Type> {
public:
  Id_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return x; }
  virtual Value_Type operator[](const Value_Type &y) { return y; }
};// end of class Id_Scaling  
template <class Value_Type>
Id_Scaling<Value_Type>::Id_Scaling(const std::string &parameter) 
{
  this->m_name="Id"; 
}

template <class Value_Type>
class Log_Scaling: public Scaling_Base<Value_Type> {
public:
  Log_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return log(x); }
  virtual Value_Type operator[](const Value_Type &y) { return exp(y); }
};// end of class Log_Scaling   
template <class Value_Type>
Log_Scaling<Value_Type>::Log_Scaling(const std::string &parameter) 
{
  this->m_name="Log"; 
}

template <class Value_Type>
class Exp_Scaling: public Scaling_Base<Value_Type> {
public:
  Exp_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return exp(x); }
  virtual Value_Type operator[](const Value_Type &y) { return log(y); }
};// end of class Exp_Scaling   
template <class Value_Type>
Exp_Scaling<Value_Type>::Exp_Scaling(const std::string &parameter) 
{
  this->m_name="Exp"; 
}

template <class Value_Type>
class Sqr_Scaling: public Scaling_Base<Value_Type> {
public:
  Sqr_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return x*x; }
  virtual Value_Type operator[](const Value_Type &y) { return sqrt(y); }
};// end of class Sqr_Scaling
template <class Value_Type>
Sqr_Scaling<Value_Type>::Sqr_Scaling(const std::string &parameter) 
{
  this->m_name="Sqr"; 
}

template <class Value_Type>
class Sqrt_Scaling: public Scaling_Base<Value_Type> {
public:
  Sqrt_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return sqrt(x); }
  virtual Value_Type operator[](const Value_Type &y) { return y*y; }
};// end of class Sqr_Scaling
template <class Value_Type>
Sqrt_Scaling<Value_Type>::Sqrt_Scaling(const std::string &parameter) 
{
  this->m_name="Sqrt"; 
}

template <class Value_Type>
class Log_B_Scaling: public Scaling_Base<Value_Type> {
private:
  Value_Type m_b, m_logb;
public:
  Log_B_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return log(x)/m_logb; }
  virtual Value_Type operator[](const Value_Type &y) { return pow(m_b,y); }
};// end of class Log_B_Scaling
template <class Value_Type>
Log_B_Scaling<Value_Type>::Log_B_Scaling(const std::string &parameter)
{
  Data_Reader reader;
  reader.SetExactMatch(false);
  reader.SetString(parameter);
  reader.ReadFromString(m_b,"Log_B_");
  m_logb=log(m_b); 
  this->m_name="Log_B_"+ToString(m_b);
}

template <class Value_Type>
class B_To_X_Scaling: public Scaling_Base<Value_Type> {
private:
  Value_Type m_b;
public:
  B_To_X_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return pow(m_b,x); }
  virtual Value_Type operator[](const Value_Type &y) 
  { return log(y)/log(m_b); }
};// end of class B_To_X_Scaling
template <class Value_Type>
B_To_X_Scaling<Value_Type>::B_To_X_Scaling(const std::string &parameter)
{
  Data_Reader reader;
  reader.SetExactMatch(false);
  reader.SetString(parameter);
  reader.ReadFromString(m_b,"B_To_X_");
  this->m_name="B_To_X_"+ToString(m_b);
}

template <class Value_Type>
class X_To_P_Scaling: public Scaling_Base<Value_Type> {
private:
  Value_Type m_p;
public:
  X_To_P_Scaling(const std::string &parameter);
  virtual Value_Type operator()(const Value_Type &x) { return pow(x,m_p); }
  virtual Value_Type operator[](const Value_Type &y) 
  { return pow(y,(Value_Type)1.0/m_p); }
};// end of class X_To_P_Scaling
template <class Value_Type>
X_To_P_Scaling<Value_Type>::X_To_P_Scaling(const std::string &parameter)
{
  Data_Reader reader;
  reader.SetExactMatch(false);
  reader.SetString(parameter);
  reader.ReadFromString(m_p,"X_To_P_");
  this->m_name="X_To_P_"+ToString(m_p);
}

template class Scaling_Base<double>;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Scaling_Base<double>
#define PARAMETER_TYPE std::string
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

template <class Class>
Scaling_Base<double> *GetVariable(const std::string &parameter) 
{
  return new Class(parameter);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Scaling_Base<double> *						\
  NAME::operator()(const std::string &parameter) const			\
  { return GetVariable<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<PRINT; }

#define DEFINE_SCALING_GETTER(CLASS,NAME,TAG,PRINT)			\
  template class CLASS;							\
  DECLARE_GETTER(NAME,TAG,Scaling_Base<double>,std::string);		\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME,PRINT) 

DECLARE_GETTER(Id_Scaling_Default_Getter,"",
	       Scaling_Base<double>,std::string);			
DEFINE_GETTER_METHOD(Id_Scaling<double>,Id_Scaling_Default_Getter)
DEFINE_PRINT_METHOD(Id_Scaling_Default_Getter,"")

DEFINE_SCALING_GETTER(Id_Scaling<double>,Id_Scaling_Getter,
		      "Id","identical")

DEFINE_SCALING_GETTER(Log_Scaling<double>,Log_Scaling_Getter,
		      "Log","logarithmical")

DEFINE_SCALING_GETTER(Exp_Scaling<double>,Exp_Scaling_Getter,
		      "Exp","exponential")

DEFINE_SCALING_GETTER(Sqr_Scaling<double>,Sqr_Scaling_Getter,
		      "Sqr","square")

DEFINE_SCALING_GETTER(Sqrt_Scaling<double>,Sqrt_Scaling_Getter,
		      "Sqrt","square root")

DEFINE_SCALING_GETTER(Log_B_Scaling<double>,Log_B_Scaling_Getter,
		      "Log_B_","logarithmical")

DEFINE_SCALING_GETTER(B_To_X_Scaling<double>,B_To_X_Scaling_Getter,
		      "B_To_X_","exponential")

DEFINE_SCALING_GETTER(X_To_P_Scaling<double>,X_To_P_Scaling_Getter,
		      "X_To_P_","power")

