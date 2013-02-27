#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Variable_Core_Scale: public Core_Scale_Setter,
			     public ATOOLS::Tag_Replacer {
  protected:

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    ATOOLS::Cluster_Amplitude *p_ampl;

  public:

    Variable_Core_Scale(const Core_Scale_Arguments &args);

    ~Variable_Core_Scale();

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc,
		  const size_t &n);

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

Variable_Core_Scale::Variable_Core_Scale
(const Core_Scale_Arguments &args): Core_Scale_Setter(args)
{
  size_t n(args.p_proc->NIn()+args.p_proc->NOut());
  p_ampl=Cluster_Amplitude::New();
  for (size_t i(0);i<n;++i)
    p_ampl->CreateLeg(Vec4D(),Flavour(kf_jet),ColorID());
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
    m_calcs.back()->AddFunction(MODEL::as->GetAIGMeanFunction());
    m_calcs.back()->SetTagReplacer(this);
    SetScale(ctag,*m_calcs.back(),n);
  }
  p_ampl->Delete();
}

Variable_Core_Scale::~Variable_Core_Scale()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
}

PDF::CParam Variable_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  p_ampl=ampl;
  double muf2(m_calcs[0]->Calculate()->Get<double>()), mur2(muf2), q2(muf2);
  if (m_calcs.size()>1) mur2=m_calcs[1]->Calculate()->Get<double>();
  if (m_calcs.size()>2) q2=m_calcs[2]->Calculate()->Get<double>();
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(q2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,q2,0.0,mur2,-1);
}

void Variable_Core_Scale::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc,const size_t &n)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): Core scale '"<<mu2tag<<"' {\n";
  msg_Indent();
  for (size_t i=0;i<n;++i) 
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

std::string Variable_Core_Scale::ReplaceTags(std::string &expr) const
{
  return m_calcs.front()->ReplaceTags(expr);
}

Term *Variable_Core_Scale::ReplaceTags(Term *term) const
{
  term->Set(p_ampl->Leg(term->Id()-100)->Mom());
  return term;
}

void Variable_Core_Scale::AssignId(Term *term)
{
  term->SetId(100+ToType<int>
	      (term->Tag().substr
	       (2,term->Tag().length()-3)));
}

DECLARE_ND_GETTER(Variable_Core_Scale,"VAR",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,Variable_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new Variable_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    Variable_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable core scale"; 
}
