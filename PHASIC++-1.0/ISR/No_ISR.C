#include "No_ISR.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace ISR;

No_ISR_at_all::No_ISR_at_all(Flavour beam)
{
  type   = std::string("(None)");
  weight = 1.;
  msg.Tracking()<<"Initialised Initial State Radiation (None) for beam "<<beam<<std::endl;
}

bool No_ISR_at_all::CalculateWeight(double x,double q2) { return 1; }
double No_ISR_at_all::Weight(Flavour fl)                { return weight; }





