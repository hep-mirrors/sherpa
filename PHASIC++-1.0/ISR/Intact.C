#include "Intact.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace ISR;

Intact::Intact(Flavour _bunch) 
{
  m_bunch  = _bunch;
  m_type   = std::string("(None)");
  m_weight = 1.;
  msg.Tracking()<<"Initialised Intact for beam "<<m_bunch<<std::endl;
}

bool Intact::CalculateWeight(double x,double q2) { return 1; }
double Intact::Weight(Flavour fl)                { return m_weight; }





