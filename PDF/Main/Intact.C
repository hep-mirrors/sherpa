#include "PDF/Main/Intact.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;

Intact::Intact(ATOOLS::Flavour _bunch):
  ISR_Base(NULL)
{
  m_bunch  = _bunch;
  m_type   = std::string("(None)");
  m_weight = 1.;
}

bool Intact::CalculateWeight(double x,double z,double kp2,double q2) 
{ 
  return 1; 
}

double Intact::Weight(ATOOLS::Flavour fl)                
{ 
  if (m_bunch.Includes(fl)) return m_weight;
  return 0.;
}





