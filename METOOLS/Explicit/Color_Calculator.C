#include "METOOLS/Explicit/Color_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Color_Calculator
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include <cassert>

using namespace METOOLS;

size_t Color_Calculator::s_cimin(1);
size_t Color_Calculator::s_cimax(3);

Color_Calculator::~Color_Calculator() {}

bool Color_Calculator::Evaluate(const CObject_Vector &j)
{
  THROW(fatal_error,"Pure virtual method called");
  return false;
}

void Color_Calculator::AddJ(CObject *const j)
{
  CObject *cc(j), *tmpcopy(NULL);

  typedef std::vector<CInfo>::const_iterator It_Type;
  for (It_Type c(m_c.begin());c!=m_c.end();++c) {

    // make a temporary copy for the following iteration (if any)
    if (c<--m_c.end()) tmpcopy=cc->Copy();
    assert(tmpcopy != NULL);

    // add cc adjusted according to c
    (*cc)(0)=c->m_cr;
    (*cc)(1)=c->m_ca;
    if (c->m_s!=Complex(1.0)) {
      if (c->m_s==Complex(-1.0)) cc->Invert();
      else if (c->m_s.imag()==0.0) cc->Multiply(c->m_s.real());
      else cc->Multiply(c->m_s);
    }
    p_v->AddJ(cc);

    // revert cc using our temporary copy
    cc=tmpcopy;
  }
}

std::ostream &METOOLS::operator<<
(std::ostream &s,const Color_Calculator::CInfo &i)
{
  return s<<"("<<i.m_cr<<","<<i.m_ca<<"|"<<i.m_s<<")";
}
