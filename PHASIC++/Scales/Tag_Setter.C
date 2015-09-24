#include "PHASIC++/Scales/Tag_Setter.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

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
  case 3:
    term->Set(p_setter->Scale(stp::res));
    return term;
  case 4:
    term->Set(sqr(p_setter->HTM()));
    return term;
  case 5:
    term->Set(sqr(p_setter->HT()));
    return term;
  case 6:
    term->Set(p_setter->PSum());
    return term;
  case 8:
    term->Set(sqr(p_setter->BeamThrust()));
    return term;
  }
  return term;
}

void Tag_Setter::AssignId(Term *term)
{
  if (term->Tag()=="MU_F2") term->SetId(1);
  else if (term->Tag()=="MU_R2") term->SetId(2);
  else if (term->Tag()=="MU_Q2") term->SetId(3);
  else if (term->Tag()=="H_TM2") term->SetId(4);
  else if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag()=="P_SUM") term->SetId(6);
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (2,term->Tag().length()-3)));
  }
}

namespace PHASIC {

  class H_TY2: public Function {
  private:

    Scale_Setter_Base *p_setter;

  public:

    inline H_TY2(Scale_Setter_Base *const setter):
      Function("H_TY2"), p_setter(setter) {}

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      double htyfac(args[0]->Get<double>()), htyexp(args[1]->Get<double>());
      Vec4D psum(0.,0.,0.,0.);
      const Vec4D_Vector &p(p_setter->Momenta());
      for (size_t i(p_setter->NIn());i<p.size();++i) psum+=p[i];
      double yboost((psum/(double)(p.size()-p_setter->NIn())).Y());
      double hty(0.0);
      for (size_t i(p_setter->NIn());i<p.size();++i) 
	hty+=p[i].PPerp()*exp(htyfac*pow(abs(p[i].Y()-yboost),htyexp));
      Term *res(Term::New(sqr(hty)));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class H_TY2

}

void Tag_Setter::SetTags(Algebra_Interpreter *const calc)
{
  calc->AddTag("MU_F2","1.0");
  calc->AddTag("MU_R2","1.0");
  calc->AddTag("MU_Q2","1.0");
  calc->AddTag("H_TM2","1.0");
  calc->AddTag("H_T2","1.0");
  calc->AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  calc->AddFunction(new H_TY2(p_setter));
  calc->AddTag("TAU_B2","1.0");
  for (size_t i=0;i<p_setter->Scales().size();++i) 
    calc->AddTag("MU_"+ToString(i)+"2","1.0");
  for (size_t i=0;i<p_setter->NIn()+p_setter->NOut();++i) 
    calc->AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
}
