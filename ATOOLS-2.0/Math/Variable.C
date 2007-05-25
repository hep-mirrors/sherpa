#include "Variable.H"

#include "Message.H"
#include "Vector.H"

using namespace ATOOLS;
  
template <class ValueType>
Variable_Base<ValueType>::Variable_Base(const std::string &name,
					const std::string &idname):
  m_name(name), m_idname(idname),
  m_selectorid(-1) 
{
  if (m_idname=="") m_idname=m_name;
}

template <class ValueType>
Variable_Base<ValueType>::~Variable_Base() {}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value
(const Vec3D *vectors,const int &n) const
{
  msg.Error()<<"Variable_Base::Value("<<vectors<<","<<n<<"): "
	     <<"Virtual method called."<<std::endl;
  return 0.0;
}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value
(const Vec4D *vectors,const int &n) const
{
  msg.Error()<<"Variable_Base::Value("<<vectors<<","<<n<<"): "
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
const std::string &Variable_Base<ValueType>::Name() const 
{
  return m_name; 
}

template <class ValueType>
const std::string &Variable_Base<ValueType>::IDName() const 
{
  return m_idname; 
}

template <class ValueType>
int Variable_Base<ValueType>::SelectorID() const 
{
  return m_selectorid; 
}

template <class ValueType>
ValueType Variable_Base<ValueType>::operator()
  (const Vec3D *vectors,const int &n) const
{ 
  return Value(vectors,n); 
}

template <class ValueType>
ValueType Variable_Base<ValueType>::operator()
  (const Vec4D *vectors,const int &n) const 
{
  return Value(vectors,n); 
}

template <class ValueType>
class No_Variable: public Variable_Base<ValueType> {
public:
  No_Variable();
};// end of class No_Variable
template <class ValueType>
No_Variable<ValueType>::No_Variable(): Variable_Base<ValueType>("") {}
  
template <class ValueType>
class Count: public Variable_Base<ValueType> {
public:
  Count();
  ValueType Value(const Vec3D *vectors,const int &n) const { return n; }
  ValueType Value(const Vec4D *vectors,const int &n) const { return n; }
};// end of class Count
template <class ValueType>
Count<ValueType>::Count(): Variable_Base<ValueType>("Count") {}
  
template <class ValueType>
class PPerp: public Variable_Base<ValueType> {
public:
  PPerp();
  ValueType Value(const Vec3D *vectors,const int &n) const
  { 
    Vec4D mom(0.0,vectors[0]);
    for (int i(1);i<n;++i) mom+=Vec4D(0.0,vectors[i]);
    return mom.PPerp(); 
  }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  {
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.PPerp();
  }
};// end of class PPerp
template <class ValueType>
PPerp<ValueType>::PPerp(): Variable_Base<ValueType>("p_\\perp","PT") 
{
  this->m_selectorid=12; 
}
  
template <class ValueType>
class EPerp: public Variable_Base<ValueType> {
public:
  EPerp();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.EPerp();
  }
};// end of class EPerp
template <class ValueType>
EPerp<ValueType>::EPerp(): Variable_Base<ValueType>("E_\\perp","ET") 
{
  this->m_selectorid=15; 
}
  
template <class ValueType>
class Energy: public Variable_Base<ValueType> {
public:
  Energy();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    double E(vectors[0][0]);
    for (int i(1);i<n;++i) E+=vectors[i][0];
    return E;
  }
};// end of class Energy
template <class ValueType>
Energy<ValueType>::Energy(): Variable_Base<ValueType>("E") 
{
  this->m_selectorid=11; 
}
  
template <class ValueType>
class Mass: public Variable_Base<ValueType> {
public:
  Mass();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.Mass();
  }
};// end of class Mass
template <class ValueType>
Mass<ValueType>::Mass(): Variable_Base<ValueType>("m") 
{
  this->m_selectorid=21; 
}
  
template <class ValueType>
class Rapidity: public Variable_Base<ValueType> {
public:
  Rapidity();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.Y();
  }
};// end of class Rapidity
template <class ValueType>
Rapidity<ValueType>::Rapidity(): Variable_Base<ValueType>("y") 
{
  this->m_selectorid=13; 
}
  
template <class ValueType>
class Eta: public Variable_Base<ValueType> {
public:
  Eta();
  ValueType Value(const Vec3D *vectors,const int &n) const 
  { 
    Vec4D mom(0.0,vectors[0]);
    for (int i(1);i<n;++i) mom+=Vec4D(0.0,vectors[i]);
    return mom.Eta(); 
  }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.Eta();
  }
};// end of class Eta
template <class ValueType>
Eta<ValueType>::Eta(): Variable_Base<ValueType>("\\eta","Eta") 
{
  this->m_selectorid=16; 
}
  
template <class ValueType>
class Theta: public Variable_Base<ValueType> {
public:
  Theta();
  ValueType Value(const Vec3D *vectors,const int &n) const 
  { 
    Vec4D mom(0.0,vectors[0]);
    for (int i(1);i<n;++i) mom+=Vec4D(0.0,vectors[i]);
    return mom.Theta(); 
  }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.Theta();
  }
};// end of class Theta
template <class ValueType>
Theta<ValueType>::Theta(): Variable_Base<ValueType>("\\theta","Theta") 
{
  this->m_selectorid=14; 
}
  
template <class ValueType>
class Phi: public Variable_Base<ValueType> {
public:
  Phi();
  ValueType Value(const Vec3D *vectors,const int &n) const 
  { 
    Vec4D mom(0.0,vectors[0]);
    for (int i(1);i<n;++i) mom+=Vec4D(0.0,vectors[i]);
    return mom.Phi(); 
  }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.Phi();
  }
};// end of class Phi
template <class ValueType>
Phi<ValueType>::Phi(): Variable_Base<ValueType>("\\phi","Phi") {}
  
template <class ValueType>
class DEta: public Variable_Base<ValueType> {
public:
  DEta();
  ValueType Value(const Vec3D *vectors) const 
  { return Vec4D(0.0,vectors[1]).DEta(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { return vectors[1].DEta(vectors[0]); }
};// end of class DEta
template <class ValueType>
DEta<ValueType>::DEta(): Variable_Base<ValueType>
("\\Delta\\eta_{ij}","DEta") {}

template <class ValueType>
class DPhi: public Variable_Base<ValueType> {
public:
  DPhi();
  ValueType Value(const Vec3D *vectors) const 
  { return Vec4D(0.0,vectors[1]).DPhi(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { return vectors[1].DPhi(vectors[0]); }
};// end of class DPhi
template <class ValueType>
DPhi<ValueType>::DPhi(): Variable_Base<ValueType>
("\\Delta\\phi_{ij}","DPhi") {}

template <class ValueType>
class DR: public Variable_Base<ValueType> {
public:
  DR();
  ValueType Value(const Vec3D *vectors) const 
  { return Vec4D(0.0,vectors[1]).DR(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { return vectors[1].DR(vectors[0]); }
};// end of class DR
template <class ValueType>
DR<ValueType>::DR(): Variable_Base<ValueType>
("\\Delta R_{ij}","DR") {}

template <class ValueType>
class Theta2: public Variable_Base<ValueType> {
public:
  Theta2();
  ValueType Value(const Vec3D *vectors) const 
  { return Vec4D(0.0,vectors[1]).Theta(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { return vectors[1].Theta(vectors[0]); }
};// end of class Theta2
template <class ValueType>
Theta2<ValueType>::Theta2(): Variable_Base<ValueType>("\\theta_{ij}","Theta2") 
{
  this->m_selectorid=22; 
}

template <class ValueType>
class T: public Variable_Base<ValueType> {
public:
  T();
};// end of class T
template <class ValueType>
T<ValueType>::T(): Variable_Base<ValueType>("t") {}
  
template <class ValueType>
class Delta: public Variable_Base<ValueType> {
public:
  Delta();
};// end of class Delta
template <class ValueType>
Delta<ValueType>::Delta(): Variable_Base<ValueType>("\\Delta(t,t_0)") {}
  
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

#define DEFINE_VARIABLE_GETTER(CLASS,NAME,TAG,PRINT)		        \
  DECLARE_GETTER(NAME,TAG,Variable_Base<double>,std::string);		\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

template class No_Variable<double>;
DEFINE_VARIABLE_GETTER(No_Variable<double>,No_Variable_Getter,"","")
template class PPerp<double>;
DEFINE_VARIABLE_GETTER(PPerp<double>,PPerp_Getter,"p_\\perp","p_\\perp")
DEFINE_VARIABLE_GETTER(PPerp<double>,PPerp_Getter_2,"PT","p_\\perp")
template class EPerp<double>;
DEFINE_VARIABLE_GETTER(EPerp<double>,EPerp_Getter,"E_\\perp","E_\\perp")
DEFINE_VARIABLE_GETTER(EPerp<double>,EPerp_Getter_2,"ET","E_\\perp")
template class Count<double>;
DEFINE_VARIABLE_GETTER(Count<double>,Count_Getter,"Count","Count")
template class Energy<double>;
DEFINE_VARIABLE_GETTER(Energy<double>,Energy_Getter,"E","E")
template class Mass<double>;
DEFINE_VARIABLE_GETTER(Mass<double>,Mass_Getter,"m","m")
template class Rapidity<double>;
DEFINE_VARIABLE_GETTER(Rapidity<double>,Rapidity_Getter,"y","y")
template class Eta<double>;
DEFINE_VARIABLE_GETTER(Eta<double>,Eta_Getter,"\\eta","\\eta")
DEFINE_VARIABLE_GETTER(Eta<double>,Eta_Getter_2,"Eta","\\eta")
template class Theta<double>;
DEFINE_VARIABLE_GETTER(Theta<double>,Theta_Getter,"\\theta","\\theta")
DEFINE_VARIABLE_GETTER(Theta<double>,Theta_Getter_2,"Theta","\\theta")
template class Phi<double>;
DEFINE_VARIABLE_GETTER(Phi<double>,Phi_Getter,"\\phi","\\phi")
DEFINE_VARIABLE_GETTER(Phi<double>,Phi_Getter_2,"Phi","\\phi")
template class Theta2<double>;
DEFINE_VARIABLE_GETTER(Theta2<double>,Theta2_Getter,
		       "\\theta_{ij}","\\theta_{ij}")
DEFINE_VARIABLE_GETTER(Theta2<double>,Theta2_Getter_2,"Theta2","\\theta_{ij}")
template class DEta<double>;
DEFINE_VARIABLE_GETTER(DEta<double>,DEta_Getter,
		       "\\Delta\\eta_{ij}","\\Delta\\eta_{ij}")
DEFINE_VARIABLE_GETTER(DEta<double>,DEta_Getter_2,"DEta","\\Delta\\eta_{ij}")
template class DPhi<double>;
DEFINE_VARIABLE_GETTER(DPhi<double>,DPhi_Getter,
		       "\\Delta\\phi_{ij}","\\Delta\\phi_{ij}")
DEFINE_VARIABLE_GETTER(DPhi<double>,DPhi_Getter_2,"DPhi","\\Delta\\phi_{ij}")
template class DR<double>;
DEFINE_VARIABLE_GETTER(DR<double>,DR_Getter,
		       "\\Delta R_{ij}","\\Delta R_{ij}")
DEFINE_VARIABLE_GETTER(DR<double>,DR_Getter_2,"DR","\\Delta R_{ij}")

template class T<double>;
DEFINE_VARIABLE_GETTER(T<double>,T_Getter,
		       "t","t")
template class Delta<double>;
DEFINE_VARIABLE_GETTER(Delta<double>,Delta_Getter,
		       "\\Delta(t,t_0)","\\Delta(t,t_0)")
