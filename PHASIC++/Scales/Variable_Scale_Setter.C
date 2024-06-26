#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace PHASIC {

  class Variable_Scale_Setter: public Scale_Setter_Base {
  protected:

    Core_Scale_Setter *p_core;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

  public:

    Variable_Scale_Setter(const Scale_Setter_Arguments &args);

    ~Variable_Scale_Setter();

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Variable_Scale_Setter,"VAR",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Variable_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Variable_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    Variable_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme";
}

Variable_Scale_Setter::Variable_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args), m_tagset(this)
{
  Settings& s = Settings::GetMainSettings();
  std::string tag(args.m_scale), core;
  size_t pos(tag.find("VAR["));
  if (pos!=std::string::npos) {
    tag=tag.substr(pos+4);
    pos=tag.find(']');
    if (pos==std::string::npos) 
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    core=tag.substr(0,pos);
    tag=tag.substr(pos+1);
  }
  if (core == "") {
    core = s["MEPS"]["CORE_SCALE"].Get<std::string>();
  }
  p_core=Core_Scale_Getter::GetObject(core,Core_Scale_Arguments(p_proc,core));
  if (p_core==NULL) THROW(fatal_error,"Invalid core scale '"+core+"'");
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
    m_calcs.back()->AddFunction(MODEL::as->GetAIGMeanFunction());
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
  delete p_core;
}

double Variable_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  if (p_proc->Flavours()[0].IsLepton()&&rpa->gen.Beam2().IsHadron()) {
    msg_Debugging()<<METHOD<<"(): Boost to Breit frame {\n";
    Vec4D pp(rpa->gen.PBeam(1)), qq(momenta[0]-momenta[2]);
    Poincare cms(pp+qq);
    double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
    Vec4D p(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
    Vec4D q(0.0,0.0,0.0,2.0*x*E);
    cms.Boost(pp);
    cms.Boost(qq);
    Poincare zrot(pp,-Vec4D::ZVEC);
    zrot.Rotate(pp);
    zrot.Rotate(qq);
    Poincare breit(p+q);
    breit.BoostBack(pp);
    breit.BoostBack(qq);
    if (!IsEqual(pp,p,1.0e-3) || !IsEqual(qq,q,1.0e-3))
      msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
    for (int i(0);i<momenta.size();++i) {
      msg_Debugging()<<"  "<<i<<": "<<m_p[i];
      cms.Boost(m_p[i]);
      zrot.Rotate(m_p[i]);
      breit.BoostBack(m_p[i]);
      msg_Debugging()<<" -> "<<m_p[i]<<"\n";
    }
    msg_Debugging()<<"}\n";
  }
  for (size_t i(0);i<m_calcs.size();++i)
    m_scale[i]=m_calcs[i]->Calculate()->Get<double>();
  for (size_t i(m_calcs.size());i<stp::size;++i) m_scale[i]=m_scale[0];
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n"
		 <<"  \\mu_q = "<<sqrt(m_scale[stp::res])<<"\n";
  for (size_t i(2);i<m_calcs.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<(p_proc?p_proc->Name():"")<<"\n";
  return m_scale[stp::fac];
}

void Variable_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<(p_proc?p_proc->Name():"")<<"' {\n";
  msg_Indent();
  m_tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}
