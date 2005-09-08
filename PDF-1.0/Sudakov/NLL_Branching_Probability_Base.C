#include "NLL_Branching_Probability_Base.H"

using namespace SHERPA;

double NLL_Branching_Probability_Base::operator()(double q, double Q) 
{
  return Gamma(q,Q);
}

double NLL_Branching_Probability_Base::operator()(double q) 
{
  return Gamma(q,m_qmax);
}

double NLL_Branching_Probability_Base::operator()()
{
  return m_defval; 
} 
