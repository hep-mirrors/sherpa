#include "PHASIC++/Scales/Tag_Setter.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

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
    term->Set(p_setter->YCut()*sqr(rpa.gen.Ecms()));
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
  else if (term->Tag()=="H_T2") term->SetId(5);
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (2,term->Tag().length()-3)));
  }
}

void Tag_Setter::SetTags(Algebra_Interpreter *const calc)
{
  calc->AddTag("MU_F2","1.0");
  calc->AddTag("MU_R2","1.0");
  calc->AddTag("Q2_CUT","1.0");
  calc->AddTag("H_T2","1.0");
  Process_Integrator *ib(p_setter->Process()->Integrator());
  ib->SetMomenta(Vec4D_Vector(ib->NIn()+ib->NOut()));
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    calc->AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
}
