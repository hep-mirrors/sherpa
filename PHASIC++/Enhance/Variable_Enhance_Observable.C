#include "PHASIC++/Enhance/Enhance_Observable_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {

  class Variable_Enhance_Observable:
    public Enhance_Observable_Base,
    public ATOOLS::Tag_Replacer {
  private:

    ATOOLS::Algebra_Interpreter m_calc;

    const ATOOLS::Vec4D *p_p;

    size_t m_n;

  public:

    Variable_Enhance_Observable(const Enhance_Arguments &args);

    double operator()(const ATOOLS::Vec4D *p,
		      const ATOOLS::Flavour *fl,
		      const size_t n);

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  };// end of class Variable_Enhance_Observable

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Variable_Enhance_Observable,"VAR",
	       Enhance_Observable_Base,Enhance_Arguments);

Enhance_Observable_Base *ATOOLS::Getter
<Enhance_Observable_Base,Enhance_Arguments,Variable_Enhance_Observable>::
operator()(const Enhance_Arguments &args) const
{
  return new Variable_Enhance_Observable(args);
}

void ATOOLS::Getter<Enhance_Observable_Base,Enhance_Arguments,Variable_Enhance_Observable>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"variable enhance observable";
}

Variable_Enhance_Observable::Variable_Enhance_Observable
(const Enhance_Arguments &args): Enhance_Observable_Base(args)
{
  std::string arg(args.m_enhance);
  size_t bpos(arg.find("VAR{")), epos(arg.find("}",bpos));
  if (bpos!=0 || epos==std::string::npos)
    THROW(fatal_error, "Invalid input in Variable_Enhance_Observable");
  arg=arg.substr(4,arg.length()-5);
  m_calc.SetTagReplacer(this);
  m_n=p_proc->NIn()+p_proc->NOut();
  p_p=&p_proc->Integrator()->Momenta().front();
  for (int i(0);i<m_n;++i)
    m_calc.AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
  m_calc.AddTag("H_T2","1.0");
  m_calc.AddTag("H_TM2","1.0");
  m_calc.AddTag("H_Tp2","1.0");
  m_calc.AddTag("H_TMp2","1.0");
  m_calc.Interprete(arg);
}

double Variable_Enhance_Observable::operator()
  (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,const size_t n)
{
  m_n=n;
  p_p=p;
  return m_calc.Calculate()->Get<double>();
}

std::string Variable_Enhance_Observable::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Variable_Enhance_Observable::ReplaceTags(Term *term) const
{
  if (term->Id()>=100) {
    if (term->Id()-100>=m_n)
      THROW(fatal_error,"p index too large");
    term->Set(p_p[term->Id()-100]);
    return term;
  }
  if (term->Id()==4) {
    double htm(0.0);
    for (size_t i(0);i<m_n;++i) htm+=p_p[i].MPerp();
    term->Set(sqr(htm));
    return term;
  }
  if (term->Id()==5) {
    double ht(0.0);
    for (size_t i(0);i<m_n;++i) ht+=p_p[i].PPerp();
    term->Set(sqr(ht));
    return term;
  }
  if (term->Id()==6) {
    double htm(0.0);
    Vec4D ewsum(0.0,0.0,0.0,0.0);
    for (size_t i(0);i<m_n;++i)
      if (!p_proc->Flavours()[i].Strong()) ewsum+=p_p[i];
      else htm+=p_p[i].MPerp();
    term->Set(sqr(htm+ewsum.MPerp()));
    return term;
  }
  if (term->Id()==7) {
    double htm(0.0);
    Vec4D ewsum(0.0,0.0,0.0,0.0);
    for (size_t i(0);i<m_n;++i)
      if (!p_proc->Flavours()[i].Strong()) ewsum+=p_p[i];
      else htm+=p_p[i].PPerp();
    term->Set(sqr(htm+ewsum.PPerp()));
    return term;
  }
  return term;
}

void Variable_Enhance_Observable::AssignId(Term *term)
{
  if (term->Tag()=="H_TM2") term->SetId(4);
  else if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag()=="H_TMp2") term->SetId(6);
  else if (term->Tag()=="H_Tp2") term->SetId(7);
  else {
    int idx(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (idx>=m_n) THROW(fatal_error,"Invalid syntax");
    term->SetId(100+idx);
  }
}
