#include "Intact.H"
#include "Message.H"

using namespace ATOOLS;
using namespace PDF;

Intact::Intact(Flavour _bunch) 
{
  m_bunch  = _bunch;
  m_type   = std::string("(None)");
  m_weight = 1.;
}

bool Intact::CalculateWeight(double x,double q2) { return 1; }
double Intact::Weight(Flavour fl)                { return m_weight; }





