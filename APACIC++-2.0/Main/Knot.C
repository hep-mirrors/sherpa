#include "Knot.H"

using std::flush;
using std::endl;

namespace APACIC {

std::ostream& operator<< (std::ostream& s, Knot& k)
{
  if (&k) {
    if (k.prev==0) s<<"#"<<flush; 
              else s<<k.prev->kn_no<<flush;
    s<<"->("<<k.kn_no<<")"<<"->("<<flush;
    if (k.left==0) s<<"#,"<<flush; 
              else s<<k.left->kn_no<<","<<flush;
    if (k.right==0) s<<"#)"<<flush; 
               else s<<k.right->kn_no<<")"<<flush;
    s<<":"<<k.part->flav()<<":"<<k.part->info()<<":"<<endl
     <<"t,tout,E2,x,z,maxpt2,thcrit :"<<k.t<<"/"<<k.tout<<"/"<<k.E2<<"/"
     <<k.x<<"/"<<k.z<<"/"<<k.maxpt2<<"/"<<k.thcrit<<", "<<endl
     <<k.part->momentum()<<","<<(k.part->momentum()).abs2()
     <<": ("<<k.part->flow(1)<<", "<<k.part->flow(2)<<") "
     <<endl;
  } 
  else { s<<"***empty knot***"<<endl; } 
  return s;
}
}
