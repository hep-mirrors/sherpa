#include "Knot.H"

using std::flush;
using std::endl;

namespace APACIC {

std::ostream& operator<< (std::ostream& s, Knot& k)
{
  if (&k) {
    if (k.prev==0) s<<"#"<<flush; 
              else s<<k.prev->kn_no<<flush;
    if (k.stat==0) s<<"->("<<k.kn_no<<")"<<"->("<<flush;
    else if (k.stat==1) s<<"->["<<k.kn_no<<"]"<<"->("<<flush;
    else s<<"->{"<<k.kn_no<<"}"<<"->("<<flush;
    if (k.left==0) s<<"#,"<<flush; 
              else s<<k.left->kn_no<<","<<flush;
    if (k.right==0) s<<"#)"<<flush; 
               else s<<k.right->kn_no<<")"<<flush;
    s<<":"<<k.part->Flav()<<":"<<k.part->Info()<<":"<<endl
     <<"t,tout,tmax,E2,x,z,maxpt2,thcrit :"<<k.t<<"/"<<k.tout<<"/"<<k.tmax<<"/"<<k.E2<<"/"
     <<k.x<<"/"<<k.z<<"/"<<k.maxpt2<<"/"<<k.thcrit<<", "<<endl
     <<k.part->Momentum()<<","<<(k.part->Momentum()).Abs2()
     <<": ("<<k.part->GetFlow(1)<<", "<<k.part->GetFlow(2)<<") "
     <<endl;
  } 
  else { s<<"***empty knot***"<<endl; } 
  return s;
}
}
