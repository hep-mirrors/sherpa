#include "No_Beam.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace BEAM;
using namespace std;

No_Beam_at_all::No_Beam_at_all(Flavour beam,double)
{
  type   = string("(None)");
  weight = 1.;
  msg.Tracking()<<"Initialised Beamhandling (None) for beam "<<beam<<endl;
}

bool No_Beam_at_all::CalculateWeight(double x,double q2) { return 1; }
double No_Beam_at_all::Weight(Flavour fl)                { return weight; }





