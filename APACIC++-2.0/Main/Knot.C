#include "Knot.H"

using std::flush;
using std::endl;

namespace APACIC {

std::ostream& operator<< (std::ostream& s, Knot& k)
{
  if (&k) {
    if (k.prev==0) s<<"#"; 
              else s<<k.prev->kn_no;
    if (k.stat==0) s<<"->("<<k.kn_no<<","<<k.didkin<<")"<<"->(";
    else if (k.stat==1) s<<"->["<<k.kn_no<<","<<k.didkin<<"]"<<"->(";
    else s<<"->{"<<k.kn_no<<","<<k.didkin<<"}"<<"->(";
    if (k.left==0) s<<"#,"; 
              else s<<k.left->kn_no<<",";
    if (k.right==0) s<<"#)"; 
               else s<<k.right->kn_no<<")";
    s<<":"<<k.part->Flav()<<":"<<k.part->Info()<<":"<<(&k)<<endl
     <<"t="<<k.t<<",tout="<<k.tout<<",tmax="<<k.tmax<<",E2="<<k.E2<<",x="
     <<k.x<<",z="<<k.z<<",maxpt2="<<k.maxpt2<<",thcrit="<<k.thcrit<<", "<<endl
     <<k.part->Momentum()<<","<<(k.part->Momentum()).Abs2()
     <<": ("<<k.part->GetFlow(1)<<", "<<k.part->GetFlow(2)<<") "
     <<endl;
    s<<" phi="<<k.phi<<" pol="<<k.polinfo<<"\n";
  } 
  else { s<<"***empty knot***"<<endl; } 
  return s;
}
}
