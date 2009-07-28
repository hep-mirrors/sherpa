#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Variable_Scale_Setter: public Scale_Setter_Base {
  protected:

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

  public:

    Variable_Scale_Setter(Process_Base *const proc,
			  const std::string &scale);

    ~Variable_Scale_Setter();

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);

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
  return new Variable_Scale_Setter(args.p_proc,args.m_scale);
}

void Variable_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme\n";
}

Variable_Scale_Setter::Variable_Scale_Setter
(Process_Base *const proc,const std::string &scale): 
  Scale_Setter_Base(proc), m_tagset(this)
{
  std::string tag(scale);
  while (true) {
    size_t pos(tag.find('{'));
    if (pos==std::string::npos) {
      if (!m_calcs.empty()) break;
      else { THROW(fatal_error,"Invalid scale '"+scale+"'"); }
    }
    tag=tag.substr(pos+1);
    pos=tag.find('}');
    if (pos==std::string::npos) THROW(fatal_error,"Invalid scale '"+scale+"'");
    std::string ctag(tag.substr(0,pos));
    tag=tag.substr(pos+1);
    m_calcs.push_back(new Algebra_Interpreter());
    m_calcs.back()->SetTagReplacer(&m_tagset);
    if (m_calcs.size()==1) m_tagset.SetCalculator(m_calcs.back());
    SetScale(ctag,*m_calcs.back());
  }
}

Variable_Scale_Setter::~Variable_Scale_Setter()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
}

double Variable_Scale_Setter::CalculateScale(const std::vector<ATOOLS::Vec4D> &momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),Max(2,(int)m_calcs.size()),0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  for (size_t i(0);i<m_calcs.size();++i)
    m_kfkey[i]=m_calcs[i]->Calculate()->Get<double>();
  m_scale[stp::ren]=m_scale[stp::fac]=m_kfkey[0];
  if (m_calcs.size()==1) m_kfkey[1]=m_kfkey[0];
  else m_scale[stp::ren]=m_kfkey[1];
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n";
  for (size_t i(0);i<m_calcs.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_kfkey[i])<<"\n";
  msg_Debugging()<<"}\n";
  return m_scale[stp::fac];
}

void Variable_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  mu2calc.AddTag("MU_F2","1.0");
  mu2calc.AddTag("MU_R2","1.0");
  mu2calc.AddTag("H_T2","1.0");
  mu2calc.AddTag("Q2_CUT","1.0");
  mu2calc.AddTag("Q2_MIN","1.0");
  Process_Integrator *ib(p_proc->Integrator());
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}

