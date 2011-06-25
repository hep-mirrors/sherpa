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
  if (term->Id()>=10) {
    if (term->Id()>=100) {
      term->Set(p_setter->Momenta()[term->Id()-100]);
      return term;
    }
    term->Set(p_setter->Scales()[term->Id()-10]);
    return term;
  }
  switch (term->Id()) {
  case 1:
    term->Set(p_setter->Scale(stp::fac));
    return term;
  case 2:
    term->Set(p_setter->Scale(stp::ren));
    return term;
  case 5:
    term->Set(sqr(p_setter->HT()));
    return term;
  case 6:
    term->Set(p_setter->PSum());
    return term;
  }
  return term;
}

void Tag_Setter::AssignId(Term *term)
{
  if (term->Tag()=="MU_F2") term->SetId(1);
  else if (term->Tag()=="MU_R2") term->SetId(2);
  else if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag()=="P_SUM") term->SetId(6);
  else if (term->Tag().find("MU_")==0) {
    term->SetId(10+ToType<int>
		(term->Tag().substr
		 (3,term->Tag().length()-4)));
  }
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
  calc->AddTag("H_T2","1.0");
  calc->AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  for (size_t i=0;i<p_setter->Scales().size();++i) 
    calc->AddTag("MU_"+ToString(i)+"2","1.0");
  Process_Integrator *ib(p_setter->Process()->Integrator());
  ib->SetMomenta(Vec4D_Vector(ib->NIn()+ib->NOut()));
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    calc->AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
}
