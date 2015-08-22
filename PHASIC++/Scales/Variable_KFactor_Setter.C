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

    ATOOLS::Algebra_Interpreter *p_calc;

    std::string m_kftag;

    ATOOLS::Vec4D_Vector m_p;
    std::vector<double>  m_mu2;

    void SetKFactor(const std::string &kftag);

  public:

    Variable_KFactor_Setter(const KFactor_Setter_Arguments &args);

    ~Variable_KFactor_Setter();

    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Variable_KFactor_Setter,"VAR",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,Variable_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new Variable_KFactor_Setter(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    Variable_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Variable kfactor scheme\n";
}

Variable_KFactor_Setter::Variable_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args)
{
  size_t pos(args.m_kfac.find('{'));
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid coupling '"+args.m_kfac+"'");
  m_kftag=args.m_kfac.substr(pos+1);
  pos=m_kftag.rfind('}');
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid coupling '"+args.m_kfac+"'");
  m_kftag=m_kftag.substr(0,pos);
  p_calc = new Algebra_Interpreter();
  p_calc->AddFunction(MODEL::as->GetAIFunction());
  p_calc->AddFunction(MODEL::aqed->GetAIFunction());
  SetKFactor(m_kftag);
  if (msg_LevelIsDebugging()) p_calc->PrintEquation();
}

Variable_KFactor_Setter::~Variable_KFactor_Setter()
{
  delete p_calc;
}

double Variable_KFactor_Setter::KFactor(const int mode) 
{
  if (!m_on) return 1.0;
  m_p=p_proc->ScaleSetter()->Momenta();
  m_mu2=p_proc->ScaleSetter()->Scales();
  return m_weight=p_calc->Calculate()->Get<double>();
}

double Variable_KFactor_Setter::KFactor(const ATOOLS::NLO_subevt &sub)
{
  if (!m_on) return 1.0;
  m_p=Vec4D_Vector(sub.p_mom,&sub.p_mom[sub.m_n]);
  m_mu2=sub.m_mu2;
  return m_weight=p_calc->Calculate()->Get<double>();
}

std::string Variable_KFactor_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *Variable_KFactor_Setter::ReplaceTags(Term *term) const
{
  switch (term->Id()) {
  case 1:
    term->Set(p_proc->ScaleSetter()->Scale(stp::ren));
    return term;
  case 2:
    term->Set(p_proc->ScaleSetter()->Scale(stp::fac));
    return term;
  case 3:
    term->Set(rpa->gen.Ecms());
    return term;
  case 4:
    term->Set(sqr(rpa->gen.Ecms()));
    return term;
  case 11:
    term->Set((double)p_proc->MaxOrder(0));
    return term;
  case 12:
    term->Set((double)p_proc->MaxOrder(1));
    return term;
  default:
    if (term->Id()>=1000) {
      term->Set(m_p[term->Id()-1000]);
      return term;
    }
    term->Set(m_mu2[term->Id()-100]);
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
  else if (term->Tag()=="Order_QCD") term->SetId(11);
  else if (term->Tag()=="Order_EW") term->SetId(12);
  else if (term->Tag().find("p[")==0) {
    term->SetId(1000+ToType<int>
		(term->Tag().substr
		 (2,term->Tag().length()-3)));
  }
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (3,term->Tag().length()-4)));
  }
}

void Variable_KFactor_Setter::SetKFactor(const std::string &kftag)
{ 
  if (kftag=="" || kftag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): coupling '"<<kftag<<"' {\n";
  msg_Indent();
  p_calc->SetTagReplacer(this);
  p_calc->AddTag("MU_F2","1.0");
  p_calc->AddTag("MU_R2","1.0");
  p_calc->AddTag("E_CMS","1.0");
  p_calc->AddTag("S_TOT","1.0");
  p_calc->AddTag("Order_QCD","0.0");
  p_calc->AddTag("Order_EW","0.0");
  if (p_proc->ScaleSetter()==NULL) THROW
    (fatal_error,"Process "+p_proc->Name()+" has no scale setter");
  for (size_t i(0);i<p_proc->ScaleSetter()->Scales().size();++i)
    p_calc->AddTag("MU_"+ToString(i)+"2","1.0");
  for (size_t i(0);i<p_proc->NIn()+p_proc->NOut();++i)
    p_calc->AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
  std::string res=p_calc->Interprete(kftag);
  msg_Debugging()<<"} -> "<<res<<"\n";
}

