#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class Variable_KFactor_Setter: 
    public KFactor_Setter_Base,
    public ATOOLS::Tag_Replacer {
  private:

    KFactor_Setter_Base *p_setter;

    ATOOLS::Algebra_Interpreter *p_calc;

    std::string m_kftag;

    void SetCoupling(const std::string &kftag);

  public:

    Variable_KFactor_Setter(Process_Base *const proc,
			    const std::string &kfac,
			    const size_t &oqcdlo,const size_t &oewlo);

    ~Variable_KFactor_Setter();

    double KFactor();

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Variable_KFactor_Setter_Getter,"VAR",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *Variable_KFactor_Setter_Getter::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new Variable_KFactor_Setter
    (args.p_proc,args.m_kfac,args.m_oqcdlo,args.m_oewlo);
}

void Variable_KFactor_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Variable kfactor scheme\n";
}

Variable_KFactor_Setter::Variable_KFactor_Setter
(Process_Base *const proc,const std::string &kfac,
 const size_t &oqcdlo,const size_t &oewlo):
  KFactor_Setter_Base(proc)
{
  size_t pos(kfac.find('{'));
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid coupling '"+kfac+"'");
  m_kftag=kfac.substr(pos+1);
  pos=m_kftag.rfind('}');
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid coupling '"+kfac+"'");
  m_kftag=m_kftag.substr(0,pos);
  p_calc = new Algebra_Interpreter();
  p_calc->AddFunction(MODEL::as->GetAIFunction());
  p_calc->AddFunction(MODEL::aqed->GetAIFunction());
}

Variable_KFactor_Setter::~Variable_KFactor_Setter()
{
  delete p_calc;
}

double Variable_KFactor_Setter::KFactor() 
{
  if (!m_on) return 1.0;
  if (!m_kfkey.Assigned()) {
     std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),2,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
    SetCoupling(m_kftag);
  }
  if (m_kfkey.Weight()!=ATOOLS::UNDEFINED_WEIGHT) return m_kfkey.Weight();
  if (p_proc->OrderQCD()<0 || p_proc->OrderEW()<0) {
    THROW(fatal_error,"Couplings not set for process '"+p_proc->Name()+"'");
  }
  m_kfkey<<p_calc->Calculate()->Get<double>();
  return m_kfkey.Weight();
}

std::string Variable_KFactor_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *Variable_KFactor_Setter::ReplaceTags(Term *term) const
{
  switch (term->Id()) {
  case 1:
    term->Set(m_kfkey.Doubles()[0]);
    return term;
  case 2:
    term->Set(m_kfkey.Doubles()[1]);
    return term;
  case 3:
    term->Set(rpa.gen.Ecms());
    return term;
  case 4:
    term->Set(sqr(rpa.gen.Ecms()));
    return term;
  default:
    term->Set(m_kfkey.Doubles()[term->Id()-100]);
    return term;
  }
  return term;
}

void Variable_KFactor_Setter::AssignId(Term *term)
{
  if (term->Tag()=="MU_R2") term->SetId(1);
  else if (term->Tag()=="MU_F2") term->SetId(2);
  else if (term->Tag()=="E_CMS") term->SetId(3);
  else if (term->Tag()=="S_TOT") term->SetId(4);
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (3,term->Tag().length()-4)));
  }
}

void Variable_KFactor_Setter::SetCoupling(const std::string &kftag)
{ 
  if (kftag=="" || kftag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): coupling '"<<kftag<<"' {\n";
  msg_Indent();
  p_calc->SetTagReplacer(this);
  p_calc->AddTag("MU_F2","1.0");
  p_calc->AddTag("MU_R2","1.0");
  p_calc->AddTag("E_CMS","1.0");
  p_calc->AddTag("S_TOT","1.0");
  for (size_t i(0);i<m_kfkey.Doubles().size();++i)
    p_calc->AddTag("MU_"+ToString(i)+"2","1.0");
  std::string res=p_calc->Interprete(kftag);
  msg_Debugging()<<"} -> "<<res<<"\n";
}

