#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class NLO_Scale_Setter;

  class NLO_Tag_Setter: public ATOOLS::Tag_Replacer {
  private:

    NLO_Scale_Setter *p_setter;

    ATOOLS::Algebra_Interpreter *p_calc;

  public:
    
    // constructor
    inline NLO_Tag_Setter(NLO_Scale_Setter *const setter): 
      p_setter(setter), p_calc(NULL) {}
    
    // member functions
    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    // inline functions
    void SetCalculator(ATOOLS::Algebra_Interpreter *const calc) { p_calc=calc; }

  };// end of class NLO_Tag_Setter

  class NLO_Scale_Setter: public Scale_Setter_Base {
  protected:

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    NLO_Tag_Setter m_muf2tagset, m_mur2tagset;

    Jet_Finder *p_jf;
    bool m_singlescale;

    ATOOLS::Info_Key m_facscalekey,m_renscalekey;
  public:

    NLO_Scale_Setter(Process_Base *const proc,
			  const std::string &mur2tag,
			  const std::string &muf2tag);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);
    double KFactor();

    void SetScale(const std::string &mu2tag,NLO_Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

    inline Jet_Finder *JF() const { return p_jf; }

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(NLO_Scale_Setter_Getter,"NLO",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *NLO_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new NLO_Scale_Setter(args.p_proc,args.m_ren,args.m_fac);
}

void NLO_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme\n";
}

NLO_Scale_Setter::NLO_Scale_Setter
(Process_Base *const proc,
 const std::string &mur2tag,const std::string &muf2tag): 
  Scale_Setter_Base(proc), m_muf2tagset(this), m_mur2tagset(this),
  p_jf(NULL), m_singlescale(0)
{
  if (muf2tag==mur2tag) m_singlescale=1;
  SetScale(muf2tag,m_muf2tagset,m_muf2calc);
  SetScale(mur2tag,m_mur2tagset,m_mur2calc);
}

std::string NLO_Tag_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *NLO_Tag_Setter::ReplaceTags(Term *term) const
{
  if (term->Tag()=="MU_F") {
    term->Set(p_setter->Scale(stp::fac));
    return term;
  }
  if (term->Tag()=="MU_R") {
    term->Set(p_setter->Scale(stp::ren));
    return term;
  }
  if (term->Tag()=="Q_CUT") {
    term->Set(p_setter->JF()->Ycut()*sqr(rpa.gen.Ecms()));
    return term;
  }
  if (term->Tag()=="Q_MIN") {
    term->Set(p_setter->JF()->ActualValue()*sqr(rpa.gen.Ecms()));
    return term;
  }

  if (term->Tag()=="H_T") {
    double ht(0.0);
    for (size_t i(p_setter->Process()->NIn());i<p_setter->Process()->ActiveMom().size();++i) 
      ht+=p_setter->Process()->ActiveMom()[i].PPerp();
    term->Set(sqr(ht));
    return term;
  }
  size_t i(atoi(term->Tag().substr(2,term->Tag().length()-3).c_str()));
  if (i>p_setter->Process()->ActiveMom().size()) THROW(fatal_error,"Momentum index too large");
  term->Set(p_setter->Process()->ActiveMom()[i]);
  return term;
}

double NLO_Scale_Setter::CalculateScale(const std::vector<ATOOLS::Vec4D> &momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfname(p_proc->Name());
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<METHOD<<"(): Assign '"<<p_proc->Name()
		   <<"' to '"<<kfname<<"','"<<kfinfo<<"'\n";
    m_kfkey.Assign(kfname,2,0,p_proc->Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);

    std::string siname=kfname.substr(kfname.rfind("__")+2);
    size_t i=siname.find("_");
    if (i!=std::string::npos) {
      siname=siname.substr(i);
      siname=kfname.substr(0,kfname.find("__"))+siname;
    }
    else siname=kfname.substr(0,kfname.find("__"));
    m_facscalekey.Assign(siname,0,0,p_proc->Integrator()->PSHandler()->GetInfo());
    m_renscalekey.Assign(siname,0,0,p_proc->Integrator()->PSHandler()->GetInfo());
  }
  if (m_facscalekey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    m_facscalekey<<m_muf2calc.Calculate()->Get<double>();
  }
  if (m_renscalekey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_singlescale) m_renscalekey<<m_facscalekey.Weight();
    else m_renscalekey<<m_mur2calc.Calculate()->Get<double>();
  }
  m_scale[stp::fac]=m_facscalekey.Weight();
  m_scale[stp::ren]=m_renscalekey.Weight();

  msg_Debugging()<<METHOD<<"(): Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<".\n";
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[1]=m_scale[stp::fac];
  return m_scale[stp::fac];
}

void NLO_Scale_Setter::SetScale
(const std::string &mu2tag,NLO_Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
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
  mu2calc.AddTag("Q_CUT","1.0");
  mu2calc.AddTag("Q_MIN","1.0");
  if ((mu2tag.find("Q_CUT")!=std::string::npos ||
       mu2tag.find("Q_MIN")!=std::string::npos) && p_jf==NULL) {
    if (p_proc->Selector()->Name()!="Combined_Selector")
      THROW(critical_error,"Q_MIN/Q_CUT implies JetFinder in selector file.");
    p_jf=(Jet_Finder*)
      ((Combined_Selector*)p_proc->Selector())
      ->GetSelector("Jetfinder");
  }
  
  for (size_t i=0;i<p_proc->ActiveMom().size();++i) {
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(p_proc->ActiveMom()[i]));
  }
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}

