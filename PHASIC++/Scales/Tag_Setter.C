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
  switch (term->Id()) {
  case 1:
    term->Set(p_setter->Scale(stp::fac));
    return term;
  case 2:
    term->Set(p_setter->Scale(stp::ren));
    return term;
  case 3:
    term->Set(JF()->Ycut()*sqr(rpa.gen.Ecms()));
    return term;
  case 4:
    term->Set(JF()->ActualValue()*sqr(rpa.gen.Ecms()));
    return term;
  case 5:
    term->Set(sqr(p_setter->HT()));
    return term;
  default:
    term->Set(p_setter->Momentum(term->Id()-100));
    return term;
  }
  return term;
}

void Tag_Setter::AssignId(Term *term)
{
  if (term->Tag()=="MU_F2") term->SetId(1);
  else if (term->Tag()=="MU_R2") term->SetId(2);
  else if (term->Tag()=="Q2_CUT") term->SetId(3);
  else if (term->Tag()=="Q2_MIN") term->SetId(4);
  else if (term->Tag()=="H_T2") term->SetId(5);
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (2,term->Tag().length()-3)));
  }
}
