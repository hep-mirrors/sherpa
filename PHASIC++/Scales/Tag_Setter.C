#include "PHASIC++/Scales/Tag_Setter.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

Jet_Finder *Tag_Setter::JF() const
{
  if (p_jf!=NULL) return p_jf;
  p_jf=(Jet_Finder*)p_setter->Process()->
    Selector()->GetSelector("Jetfinder");
  if (p_jf==NULL) THROW(critical_error,"Jet finder not found");
  return p_jf;
}

std::string Tag_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *Tag_Setter::ReplaceTags(Term *term) const
{
  if (term->Tag()=="MU_F") {
    term->Set(p_setter->Scale(stp::fac));
    return term;
  }
  if (term->Tag()=="MU_R") {
    term->Set(p_setter->Scale(stp::ren));
    return term;
  }
  if (term->Tag()=="Q2_CUT") {
    term->Set(JF()->Ycut()*sqr(rpa.gen.Ecms()));
    return term;
  }
  if (term->Tag()=="Q2_MIN") {
    term->Set(JF()->ActualValue()*sqr(rpa.gen.Ecms()));
    return term;
  }
  if (term->Tag()=="H_T") {
    term->Set(sqr(p_setter->HT()));
    return term;
  }
  size_t i(ToType<size_t>(term->Tag().substr(2,term->Tag().length()-3)));
  term->Set(p_setter->Momentum(i));
  return term;
}

