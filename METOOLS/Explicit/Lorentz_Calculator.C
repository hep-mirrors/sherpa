#include "METOOLS/Explicit/Lorentz_Calculator.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Lorentz_Calculator
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Single_Vertex.H"

using namespace METOOLS;

Lorentz_Calculator::Lorentz_Calculator(const Vertex_Key &key): 
  p_v(key.p_v)
{
  m_r.resize(key.p_mv->in.size());
  for (size_t i(0);i<key.p_mv->in.size();++i)
    m_r[i]=key.p_mv->in[i].Majorana()?0:(key.p_mv->in[i].IsAnti()?-1:1);
}

Lorentz_Calculator::~Lorentz_Calculator() {}

void Lorentz_Calculator::Evaluate()
{
  THROW(fatal_error,"Invalid call");
}
