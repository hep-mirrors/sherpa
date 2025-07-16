#include "PDF/Main/Intact.H"

using namespace PDF;

Intact::Intact(const ATOOLS::Flavour& _bunch):
  ISR_Base(nullptr)
{
  m_bunch  = _bunch;
  m_type   = isrtype::intact;
  m_weight = 1.;
}

bool Intact::CalculateWeight(double x,double z,double kp2,double q2,int warn) 
{ 
  return true;
}

double Intact::Weight(ATOOLS::Flavour fl)                
{ 
  if (m_bunch.Includes(fl)) return m_weight;
  return 0.;
}





