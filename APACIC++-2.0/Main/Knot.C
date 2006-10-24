#include "Knot.H"

#include "Message.H"
#include "Blob.H"

using namespace ATOOLS;
using namespace APACIC;

std::ostream &APACIC::operator<<(std::ostream& s,const Knot &k)
{
  if (&k) {
    if (k.prev==0) s<<"#"; 
    else s<<k.prev->kn_no;
    if (k.stat==0) s<<"->("<<k.kn_no<<","<<(k.decay?k.decay->kn_no:0)
		    <<","<<k.shower<<")"<<"->(";
    else if (k.stat==1) s<<"->["<<k.kn_no<<","<<(k.decay?k.decay->kn_no:0)
			 <<","<<k.shower<<"]"<<"->(";
    else s<<"->{"<<k.kn_no<<","<<(k.decay?k.decay->kn_no:0)
	  <<","<<k.shower<<"}"<<"->(";
    if (k.left==0) s<<"#,"; 
    else s<<k.left->kn_no<<",";
    if (k.right==0) s<<"#)"; 
    else s<<k.right->kn_no<<")";
    s<<":"<<k.part->Flav()<<":"<<k.part->Info()<<":"<<(&k)<<"\n"
     <<"tmo="<<k.tmo<<"("<<sqrt(dabs(k.tmo))<<"),t="<<k.t<<"("<<sqrt(dabs(k.t))
     <<"),tout="<<k.tout<<"("<<sqrt(dabs(k.tout))<<"),tmax="<<k.tmax<<"("
     <<sqrt(dabs(k.tmax))<<"),E="<<sqrt(k.E2)<<",z="<<k.z<<",zs="<<k.zs
     <<",x="<<k.x<<",\nptlcm="<<sqrt(k.pt2lcm)<<",maxpt="<<sqrt(k.maxpt2)
     <<",smaxpt="<<sqrt(k.smaxpt2)
     <<",thcrit="<<k.thcrit<<",sthcrit="<<k.sthcrit<<",minpt="<<sqrt(k.minpt2)
     <<",qjv="<<k.qjv<<",qljv="<<k.qljv<<",maxjets="<<k.maxjets<<",\n"
     <<k.part->Momentum()<<","<<(k.part->Momentum()).Abs2()
     <<": ("<<k.part->GetFlow(1)<<", "<<k.part->GetFlow(2)<<") {"
     <<(k.part->ProductionBlob()?k.part->ProductionBlob()->Id():0)<<","
     <<(k.part->DecayBlob()?k.part->DecayBlob()->Id():0)<<"} "<<k.didkin<<"\n";
    s<<" phi="<<k.phi<<" pol="<<k.polinfo<<"\n";
  } 
  else { 
    s<<"***empty knot***\n"; 
  } 
  return s;
}

Knot::Knot():
  prev(NULL), left(NULL), right(NULL), decay(NULL),
  part(NULL), stat(0), kn_no(-1),
  shower(1), maxjets(1000), didkin(false), 
  t(0.0), tout(0.0), tmax(0.0), z(0.0), zs(0.0),
  E2(0.0), costh(0.0), phi(0.0), thcrit(M_PI),
  maxpt2(1.0e10), x(0.), pt2lcm(0.0), smaxpt2(1.0e10), sthcrit(M_PI),
  minpt2(0.0), qjv(1.0e10), qljv(0.0), tmo(0.0), lz(0.0), lE2(0.0) {}

Knot::Knot(Knot * k):
  prev(k->prev), left(k->left), right(k->right), decay(k->decay), 
  part(new ATOOLS::Particle(*k->part)), lp(k->lp), 
  stat(k->stat), kn_no(k->kn_no),
  shower(k->shower), maxjets(k->maxjets), didkin(k->didkin), 
  t(k->t), tout(k->tout), tmax(k->tmax), z(k->z), zs(0.0),
  E2(k->E2), costh(k->costh), phi(k->phi), thcrit(k->thcrit), 
  maxpt2(k->maxpt2), x(k->x), pt2lcm(0.0), smaxpt2(k->smaxpt2), 
  sthcrit(k->sthcrit), minpt2(k->minpt2), qjv(k->qjv), qljv(k->qljv), 
  tmo(k->tmo), lz(k->lz), lE2(k->lE2), polinfo(k->polinfo) 
{
  part->SetProductionBlob(k->part->ProductionBlob());
  part->SetDecayBlob(k->part->DecayBlob());
}

void Knot::CopyData(const Knot *const k) 
{
  if (part!=NULL && k->part) {
    part->SetProductionBlob(NULL);
    part->SetDecayBlob(NULL);
    part->Copy(k->part);
  }
  lp=k->lp;
  stat=k->stat;
  kn_no=k->kn_no;
  didkin=k->didkin;
  shower=k->shower;
  maxjets=k->maxjets;
  tmo=k->tmo;
  lz=k->lz;
  lE2=k->lE2;
  t=k->t;
  tout=k->tout;
  tmax=k->tmax;
  E2=k->E2;
  z=k->z;
  zs=k->zs;
  x=k->x;
  maxpt2=k->maxpt2;
  pt2lcm=k->pt2lcm;
  costh=k->costh;
  phi=k->phi;
  thcrit=k->thcrit;
  sthcrit=k->sthcrit;
  smaxpt2=k->smaxpt2;
  minpt2=k->minpt2;
  qjv=k->qjv;
  qljv=k->qljv;
  polinfo=k->polinfo;
  decay=k->decay;
}

void Knot::Copy(const Knot *const k) 
{
  prev=k->prev;
  left=k->left;
  right=k->right;
  CopyData(k);
}

void Knot::Store()
{
  lp=part->Momentum();
  lz=z;
  lE2=E2;
  if (left) {
    left->Store();
    right->Store();
  }
}

void Knot::Restore()
{
  if (lp!=Vec4D()) {
    part->SetMomentum(lp);
    z=lz;
    E2=lE2;
  }
  if (left) {
    left->Restore();
    right->Restore();
  }
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
