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

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    Tag_Setter m_muf2tagset, m_mur2tagset;

  public:

    Variable_Scale_Setter(Process_Base *const proc,
			  const std::string &mur2tag,
			  const std::string &muf2tag);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
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
  return new Variable_Scale_Setter(args.p_proc,args.m_ren,args.m_fac);
}

void Variable_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme\n";
}

Variable_Scale_Setter::Variable_Scale_Setter
(Process_Base *const proc,
 const std::string &mur2tag,const std::string &muf2tag): 
  Scale_Setter_Base(proc), m_muf2tagset(this), m_mur2tagset(this)
{
  SetScale(muf2tag,m_muf2tagset,m_muf2calc);
  SetScale(mur2tag,m_mur2tagset,m_mur2calc);
}

double Variable_Scale_Setter::CalculateScale(const std::vector<ATOOLS::Vec4D> &momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),2,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  m_scale[stp::fac]=m_muf2calc.Calculate()->Get<double>();
  m_scale[stp::ren]=m_mur2calc.Calculate()->Get<double>();
  msg_Debugging()<<METHOD<<"(): Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<".\n";
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[1]=m_scale[stp::fac];
  return m_scale[stp::fac];
}

void Variable_Scale_Setter::SetScale
(const std::string &mu2tag,Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  mu2tagset.SetCalculator(&mu2calc);
  mu2calc.SetTagReplacer(&mu2tagset);
  mu2calc.AddTag("MU_F","1.0");
  mu2calc.AddTag("MU_R","1.0");
  mu2calc.AddTag("H_T","1.0");
  mu2calc.AddTag("Q2_CUT","1.0");
  mu2calc.AddTag("Q2_MIN","1.0");
  Process_Integrator *ib(p_proc->Integrator());
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}

