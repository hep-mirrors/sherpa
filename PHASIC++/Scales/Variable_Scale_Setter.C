#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Variable_Scale_Setter: public Scale_Setter_Base {
  protected:

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

  public:

    Variable_Scale_Setter(const Scale_Setter_Arguments &args);

    ~Variable_Scale_Setter();

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p,
			  const int mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Variable_Scale_Setter_Getter,"VAR",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Variable_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Variable_Scale_Setter(args);
}

void Variable_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme";
}

Variable_Scale_Setter::Variable_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args), m_tagset(this)
{
  std::string tag(args.m_scale);
  while (true) {
    size_t pos(tag.find('{'));
    if (pos==std::string::npos) {
      if (!m_calcs.empty()) break;
      else { THROW(fatal_error,"Invalid scale '"+args.m_scale+"'"); }
    }
    tag=tag.substr(pos+1);
    pos=tag.find('}');
    if (pos==std::string::npos) 
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    std::string ctag(tag.substr(0,pos));
    tag=tag.substr(pos+1);
    m_calcs.push_back(new Algebra_Interpreter());
    m_calcs.back()->SetTagReplacer(&m_tagset);
    if (m_calcs.size()==1) m_tagset.SetCalculator(m_calcs.back());
    SetScale(ctag,*m_calcs.back());
  }
  m_scale.resize(Max(m_scale.size(),m_calcs.size()));
  SetCouplings();
}

Variable_Scale_Setter::~Variable_Scale_Setter()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
}

double Variable_Scale_Setter::CalculateScale
(const std::vector<ATOOLS::Vec4D> &momenta,const int mode) 
{
  if (mode==1) return m_scale[stp::fac];
  if (m_escale.size()) {
    m_scale[stp::fac]=m_escale[stp::fac];
    m_scale[stp::ren]=m_escale[stp::ren];
    p_cpls->Calculate();
    return m_scale[stp::fac];    
  }
  for (size_t i(0);i<m_calcs.size();++i)
    m_scale[i]=m_calcs[i]->Calculate()->Get<double>();
  if (m_calcs.size()==1) m_scale[1]=m_scale[0];
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n";
  for (size_t i(2);i<m_calcs.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<p_proc->Name()<<"\n";
  p_cpls->Calculate();
  return m_scale[stp::fac];
}

void Variable_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  m_tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

