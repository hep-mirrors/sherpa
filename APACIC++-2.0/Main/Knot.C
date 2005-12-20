#include "Knot.H"

#include "Message.H"

using namespace ATOOLS;
using namespace APACIC;

std::ostream &APACIC::operator<<(std::ostream& s,const Knot &k)
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
    s<<":"<<k.part->Flav()<<":"<<k.part->Info()<<":"<<(&k)<<"\n"
     <<"t="<<k.t<<",tout="<<k.tout<<",tmax="<<k.tmax<<",E2="<<k.E2<<",x="
     <<k.x<<",z="<<k.z<<",zs="<<k.zs<<",maxpt2="<<k.maxpt2<<",smaxpt2="
     <<k.smaxpt2<<",pt2lcm="<<k.pt2lcm
     <<",thcrit="<<k.thcrit<<",sthcrit="<<k.sthcrit<<",\n"<<k.part->Momentum()<<","
     <<(k.part->Momentum()).Abs2()
     <<": ("<<k.part->GetFlow(1)<<", "<<k.part->GetFlow(2)<<")\n";
    s<<" phi="<<k.phi<<" pol="<<k.polinfo<<"\n";
  } 
  else { 
    s<<"***empty knot***\n"; 
  } 
  return s;
}

bool Knot::CheckMomentumConservation() const 
{
  if (left==NULL || stat==3) return true;
  bool success(true);
  Vec4D p(part->Momentum());
  Vec4D p1(left->part->Momentum()), p2(right->part->Momentum());
  if (!(p==p1+p2)) {
    msg.Error()<<METHOD<<"(): Four momentum not conserved in knot "
	       <<kn_no<<"\n   p      = "<<p<<"\n   p_miss = "<<(p-p1-p2)
	       <<"\n   p1     = "<<p1<<"\n   p2     = "<<p2<<std::endl;
    success=false;
  }
  return success;
}
