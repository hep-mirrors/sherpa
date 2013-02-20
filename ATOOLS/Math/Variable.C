#include "ATOOLS/Math/Variable.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
  
template <class ValueType>
Variable_Base<ValueType>::Variable_Base(const std::string &name,
					const std::string &idname):
  m_name(name), m_idname(idname)
{
  if (m_idname=="") m_idname=m_name;
}

template <class ValueType>
Variable_Base<ValueType>::~Variable_Base() {}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value
(const Vec3D *vectors,const int &n) const
{
  msg_Error()<<"Variable_Base::Value("<<vectors<<","<<n<<"): "
	     <<"Virtual method called."<<std::endl;
  return 0.0;
}

template <class ValueType>
ValueType Variable_Base<ValueType>::Value
(const Vec4D *vectors,const int &n) const
{
  msg_Error()<<"Variable_Base::Value("<<vectors<<","<<n<<"): "
	     <<"Virtual method called."<<std::endl;
  return 0.0;
}

template <class ValueType>
bool Variable_Base<ValueType>::Init(const std::string &name)
{
  return true;
}

template <class ValueType>
Algebra_Interpreter *Variable_Base<ValueType>::GetInterpreter() const
{
  return NULL;
}

template <class ValueType>
void Variable_Base<ValueType>::ShowVariables(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<"Variable_Base::ShowVariables(): {\n\n";
  Variable_Getter::PrintGetterInfo(msg->Out(),20);
  msg_Out()<<"\n}"<<std::endl;
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
std::string Variable_Base<ValueType>::SelectorID() const 
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

double Get(Algebra_Interpreter *const inter)
{
  Term *res(inter->Calculate());
  return res->Get<double>();
}
  
template <class ValueType>
class Calc_Variable: public Variable_Base<ValueType>,
		     public Tag_Replacer {
private:
  std::string m_formula;
  Algebra_Interpreter *p_interpreter;
  Tag_Replacer *p_replacer;
  mutable std::vector<Vec4D> m_p;
public:
  Calc_Variable(const std::string &tag);
  ~Calc_Variable();
  ValueType Value(const Vec3D *vectors,const int &n) const 
  { 
    m_p.resize(n);
    for (int i(0);i<n;++i) m_p[i]=Vec4D(0.0,vectors[i]);
    return Get(p_interpreter);
  }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    m_p.resize(n);
    for (int i(0);i<n;++i) m_p[i]=vectors[i];
    return Get(p_interpreter);
  }
  Algebra_Interpreter *GetInterpreter() const
  {
    return p_interpreter;
  }
  bool Init(const std::string &name);
  std::string ReplaceTags(std::string &expr) const;
  ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
  void AssignId(ATOOLS::Term *term);
};// end of class Calc_Variable
template <class ValueType>
Calc_Variable<ValueType>::Calc_Variable(const std::string &tag): 
  Variable_Base<ValueType>("Calc"), m_formula(tag), 
  p_interpreter(new Algebra_Interpreter()), p_replacer(NULL)
{
  p_interpreter->SetTagReplacer(this);
  Init(m_formula);
}
template <class ValueType> bool
Calc_Variable<ValueType>::Init(const std::string &name)
{
  m_formula=name;
  msg_Debugging()<<METHOD<<"(): m_formula = '"<<m_formula<<"'\n";
  size_t cbpos(m_formula.find("{"));
  if (cbpos!=std::string::npos) {
    std::string reps(m_formula.substr(cbpos+1));
    m_formula=m_formula.substr(0,cbpos);
    if ((cbpos=reps.rfind("}"))==std::string::npos) 
      THROW(fatal_error,"Invalid syntax");
    reps=reps.substr(0,cbpos);
    p_replacer=dynamic_cast<Tag_Replacer*>
      ((Tag_Replacer*)ToType<long unsigned int>(reps));
    if (p_replacer==NULL) THROW(fatal_error,"Invalid pointer");
  }
  size_t bpos(m_formula.find("("));
  if (bpos==std::string::npos) return false;
  m_formula=m_formula.substr(bpos);
  if ((bpos=m_formula.rfind(")"))==std::string::npos) return false;
  m_formula=m_formula.substr(1,bpos-1);
  if (m_formula.length()==0) return false;
  size_t pos(m_formula.find("p[")), nmf(0);
  while (pos!=std::string::npos) {
    std::string ex(m_formula.substr(pos+2,m_formula.find("]",pos)-pos-2));
    p_interpreter->AddTag("p["+ex+"]","(1.0,0.0,0.0,1.0)");
    pos=m_formula.find("p[",pos+ex.length()+1);
    nmf=Max(nmf,ToType<size_t>(ex));
  }
  m_p.resize(nmf+1);
  p_interpreter->Interprete(m_formula);
  if (msg_LevelIsTracking()) p_interpreter->PrintEquation();
  return true;
}
template <class ValueType>
Calc_Variable<ValueType>::~Calc_Variable()
{
  delete p_interpreter;
}
template <class ValueType>
std::string Calc_Variable<ValueType>::ReplaceTags(std::string &expr) const
{
  return p_interpreter->ReplaceTags(expr);
}
template <class ValueType>
ATOOLS::Term *Calc_Variable<ValueType>::ReplaceTags(ATOOLS::Term *term) const
{
  if (term->Id()>=100) term->Set(m_p[term->Id()-100]);
  else if (p_replacer!=NULL) return p_replacer->ReplaceTags(term);
  else THROW(fatal_error,"Invalid tag.");
  return term;
}
template <class ValueType>
void Calc_Variable<ValueType>::AssignId(ATOOLS::Term *term)
{
  if (term->Tag().find("p[")==0) {
    size_t i(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (i>=m_p.size()) THROW(fatal_error,"Invalid tag.");
    term->SetId(100+i);
  }
  else if (p_replacer!=NULL) p_replacer->AssignId(term);
  else THROW(fatal_error,"Invalid tag.");
}
  
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
  this->m_selectorid="PT"; 
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
  this->m_selectorid="ET"; 
}
  
template <class ValueType>
class MPerp: public Variable_Base<ValueType> {
public:
  MPerp();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    Vec4D mom(vectors[0]);
    for (int i(1);i<n;++i) mom+=vectors[i];
    return mom.MPerp();
  }
};// end of class MPerp
template <class ValueType>
MPerp<ValueType>::MPerp(): Variable_Base<ValueType>("m_\\perp","mT") {}

template <class ValueType>
class MTWW: public Variable_Base<ValueType> {
public:
  MTWW();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    if (n!=3) THROW(fatal_error,"Variable MTWW only defined for three momenta."); 
    Vec4D ll(vectors[0]+vectors[1]);
    double m2ll(ll*ll);
    Vec4D miss(vectors[2]);
    Vec3D llp = Vec3D(ll[1],ll[2],0.);
    Vec3D missp = Vec3D(miss[1],miss[2],0.);
    return sqrt(sqr(sqrt(llp*llp+m2ll)+sqrt(missp*missp+m2ll))-(llp+missp)*(llp+missp));
  }
};// end of class MTWW
template <class ValueType>
MTWW<ValueType>::MTWW(): Variable_Base<ValueType>("m_\\perp(WW)","mTWW") {}

  
template <class ValueType>
class HT: public Variable_Base<ValueType> {
public:
  HT();
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { 
    double ht(vectors[0].PPerp());
    for (int i(1);i<n;++i) ht+=vectors[i].PPerp();
    return ht;
  }
};// end of class HT
template <class ValueType>
HT<ValueType>::HT(): Variable_Base<ValueType>("H_T","HT") {}
  
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
  this->m_selectorid="Energy"; 
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
  this->m_selectorid="Mass"; 
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
  this->m_selectorid="Rapidity"; 
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
  this->m_selectorid="PseudoRapidity"; 
}
  
template <class ValueType>
class BTheta: public Variable_Base<ValueType> {
public:
  BTheta();
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
BTheta<ValueType>::BTheta(): Variable_Base<ValueType>("\\theta","Theta") 
{
  this->m_selectorid="BeamAngle"; 
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
class DY: public Variable_Base<ValueType> {
public:
  DY();
  ValueType Value(const Vec3D *vectors) const 
  { return Vec4D(0.0,vectors[1]).DY(Vec4D(0.0,vectors[0])); }
  ValueType Value(const Vec4D *vectors,const int &n) const 
  { return vectors[1].DY(vectors[0]); }
};// end of class DY
template <class ValueType>
DY<ValueType>::DY(): Variable_Base<ValueType>
("\\Delta y_{ij}","DY") {}

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
  this->m_selectorid="Angle"; 
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
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

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

#define DEFINE_VARIABLE_GETTER(CLASS,NAME,TAG,PRINT,DISP)		\
  DECLARE_ND_GETTER(NAME,TAG,Variable_Base<double>,std::string,DISP);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

template class No_Variable<double>;
DEFINE_VARIABLE_GETTER(No_Variable<double>,No_Variable_Getter,"","",0)

template class Calc_Variable<double>;
DECLARE_ND_GETTER(Calc_Variable_Getter,"Calc",
		  Variable_Base<double>,std::string,1);
Variable_Base<double> *Calc_Variable_Getter::operator()
(const std::string &parameter) const			
{ return new Calc_Variable<double>(parameter); }
void Calc_Variable_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ str<<"calculator, usage: Calc(<formula>)"; }

template class PPerp<double>;
DEFINE_VARIABLE_GETTER(PPerp<double>,PPerp_Getter,"p_\\perp","p_\\perp",0)
DEFINE_VARIABLE_GETTER(PPerp<double>,PPerp_Getter_2,"PT","p_\\perp",1)
template class EPerp<double>;
DEFINE_VARIABLE_GETTER(EPerp<double>,EPerp_Getter,"E_\\perp","E_\\perp",0)
DEFINE_VARIABLE_GETTER(EPerp<double>,EPerp_Getter_2,"ET","E_\\perp",1)
template class MPerp<double>;
DEFINE_VARIABLE_GETTER(MPerp<double>,MPerp_Getter,"m_\\perp","m_\\perp",0)
DEFINE_VARIABLE_GETTER(MPerp<double>,MPerp_Getter_2,"mT","m_\\perp",1)
template class MTWW<double>;
DEFINE_VARIABLE_GETTER(MTWW<double>,MTWW_Getter,"m_\\perp(WW)","m_\\perp(WW)",0)
DEFINE_VARIABLE_GETTER(MTWW<double>,MTWW_Getter_2,"mTWW","m_\\perp(WW)",1)
template class HT<double>;
DEFINE_VARIABLE_GETTER(HT<double>,HT_Getter_2,"H_T","H_T",0)
DEFINE_VARIABLE_GETTER(HT<double>,HT_Getter_3,"HT","H_T",1)
template class Count<double>;
DEFINE_VARIABLE_GETTER(Count<double>,Count_Getter,"Count","Count",0)
DEFINE_VARIABLE_GETTER(Count<double>,Count_Getter_2,"N","number",1)
template class Energy<double>;
DEFINE_VARIABLE_GETTER(Energy<double>,Energy_Getter,"E","E",1)
template class Mass<double>;
DEFINE_VARIABLE_GETTER(Mass<double>,Mass_Getter,"m","m",1)
template class Rapidity<double>;
DEFINE_VARIABLE_GETTER(Rapidity<double>,Rapidity_Getter,"y","y",1)
template class Eta<double>;
DEFINE_VARIABLE_GETTER(Eta<double>,Eta_Getter,"\\eta","\\eta",0)
DEFINE_VARIABLE_GETTER(Eta<double>,Eta_Getter_2,"Eta","\\eta",1)
template class BTheta<double>;
DEFINE_VARIABLE_GETTER(BTheta<double>,BTheta_Getter,"\\theta","\\theta",0)
DEFINE_VARIABLE_GETTER(BTheta<double>,BTheta_Getter_2,"Theta","\\theta",1)
template class Phi<double>;
DEFINE_VARIABLE_GETTER(Phi<double>,Phi_Getter,"\\phi","\\phi",0)
DEFINE_VARIABLE_GETTER(Phi<double>,Phi_Getter_2,"Phi","\\phi",1)
template class Theta2<double>;
DEFINE_VARIABLE_GETTER(Theta2<double>,Theta2_Getter,
		       "\\theta_{ij}","\\theta_{ij}",0)
DEFINE_VARIABLE_GETTER(Theta2<double>,Theta2_Getter_2,
		       "Theta2","\\theta_{ij}",1)
template class DEta<double>;
DEFINE_VARIABLE_GETTER(DEta<double>,DEta_Getter,
		       "\\Delta\\eta_{ij}","\\Delta\\eta_{ij}",0)
DEFINE_VARIABLE_GETTER(DEta<double>,DEta_Getter_2,
		       "DEta","\\Delta\\eta_{ij}",1)
template class DY<double>;
DEFINE_VARIABLE_GETTER(DY<double>,DY_Getter,
		       "\\Delta y_{ij}","\\Delta y_{ij}",0)
DEFINE_VARIABLE_GETTER(DY<double>,DY_Getter_2,
		       "DY","\\Delta y_{ij}",1)
template class DPhi<double>;
DEFINE_VARIABLE_GETTER(DPhi<double>,DPhi_Getter,
		       "\\Delta\\phi_{ij}","\\Delta\\phi_{ij}",0)
DEFINE_VARIABLE_GETTER(DPhi<double>,DPhi_Getter_2,
		       "DPhi","\\Delta\\phi_{ij}",1)
template class DR<double>;
DEFINE_VARIABLE_GETTER(DR<double>,DR_Getter,
		       "\\Delta R_{ij}","\\Delta R_{ij}",0)
DEFINE_VARIABLE_GETTER(DR<double>,DR_Getter_2,"DR","\\Delta R_{ij}",1)

template class T<double>;
DEFINE_VARIABLE_GETTER(T<double>,T_Getter,"t","t",0)
template class Delta<double>;
DEFINE_VARIABLE_GETTER(Delta<double>,Delta_Getter,
		       "\\Delta(t,t_0)","\\Delta(t,t_0)",0)
