#include "MCATNLO/Main/CS_Cluster_Definitions.H"

#include "MCATNLO/Showers/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace MCATNLO;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

const double s_uxeps=1.0e-3;

CS_Cluster_Definitions::CS_Cluster_Definitions
(Shower *const shower,const int kmode):
  p_shower(shower), m_mode(0), m_kmode(kmode), m_amode(0) {}

CParam CS_Cluster_Definitions::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin,const int kmode)
{
  m_mode=m_kmode;
  CS_Parameters cs(KT2(&ampl,ampl.Leg(i),ampl.Leg(j),ampl.Leg(k),mo,ms,kmode));
  m_mode=0;
  return CParam(cs.m_kt2,cs.m_ws,cs.m_x,cs.m_mu2,cs.m_kin,cs.m_kmode);
}

CS_Parameters CS_Cluster_Definitions::KT2
(const ATOOLS::Cluster_Amplitude *ampl,
 const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,
 ATOOLS::Mass_Selector *const ms,const int ikin,const int kmode)
{
  p_ms=ms;
  int kin(ikin<0?p_shower->KinScheme():ikin), col(1);
  if ((i->Id()&3)<(j->Id()&3)) {
    std::swap<const Cluster_Leg*>(i,j);
    col=-1;
  }
  p_b=ampl->Leg(i==ampl->Leg(0)?1:0);
  Vec4D pi(i->Mom()), pj(j->Mom()), pk(k->Mom());
  double Q2=(pi+pj+pk).Abs2(), mb2=p_ms->Mass2(p_b->Flav());
  double mi2=p_ms->Mass2(i->Flav()), mj2=p_ms->Mass2(j->Flav());
  double mk2=p_ms->Mass2(k->Flav()), mij2=p_ms->Mass2(mo);
  if (!(i->Id()&3) && mi2>10.0 && !i->Flav().Strong()) mi2=pi.Abs2();
  if (!(j->Id()&3) && mj2>10.0 && !j->Flav().Strong()) mj2=pj.Abs2();
  if (!(k->Id()&3) && mk2>10.0 && !k->Flav().Strong()) mk2=pk.Abs2();
  if (!(i->Id()&3) && !(j->Id()&3) && mij2>10.0 && !mo.Strong()) {
    mij2=(pi+pj).Abs2();
    pk[0]=pk[0]<0.0?-pk.PSpat():pk.PSpat();
    mk2=0.0;
  }
  CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),
		   1.0,1.0,0.0,0.0,0.0,
		   ((i->Id()&3)?1:0)|((k->Id()&3)?2:0),kin);
  cs.m_wk=cs.m_ws=-1.0;
  if ((i->Id()&3)==0) {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) {
	Kin_Args ff(ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,1|(kin?4:0)));
	if (ff.m_stat!=1) return cs;
	double kt2=2.0*(pi*pj)*ff.m_z*(1.0-ff.m_z)
	  -sqr(1.0-ff.m_z)*mi2-sqr(ff.m_z)*mj2;
	cs=CS_Parameters(kt2,ff.m_z,ff.m_y,ff.m_phi,1.0,Q2,0,kin);
      }
      else {
	Kin_Args fi(ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,1|(kin?4:0)));
	Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
	if (fi.m_pk.PPlus()>sum.PPlus() || fi.m_y>1.0 ||
	    fi.m_pk.PMinus()>sum.PMinus() || fi.m_stat!=1) return cs;
	double kt2=2.0*(pi*pj)*fi.m_z*(1.0-fi.m_z)
	  -sqr(1.0-fi.m_z)*mi2-sqr(fi.m_z)*mj2;
	cs=CS_Parameters(kt2,fi.m_z,fi.m_y,fi.m_phi,1.0-fi.m_y,Q2,2,kin);
      }
    }
  }
  else {
    if ((j->Id()&3)==0) {
      Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
      if ((k->Id()&3)==0) {
	Kin_Args fi(ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-p_b->Mom(),1|(kin?4:0)));
	if (fi.m_pi.PPlus()>sum.PPlus() || fi.m_z<0.0 ||
	    fi.m_pi.PMinus()>sum.PMinus() || fi.m_stat!=1) return cs;
	double kt2=-2.0*(pi*pj)*(1.0-fi.m_z)-mj2-sqr(1.0-fi.m_z)*mi2;
	cs=CS_Parameters(kt2,fi.m_z,fi.m_y,fi.m_phi,fi.m_z,Q2,1,fi.m_mode);
      }
      else {
	Kin_Args ii(ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,1|(kin?4:0)));
	if (ii.m_pi.PPlus()>sum.PPlus() || ii.m_z<0.0 ||
	    ii.m_pi.PMinus()>sum.PMinus() || ii.m_stat!=1) return cs;
	double kt2=-2.0*(pi*pj)*(1.0-ii.m_z)-mj2-sqr(1.0-ii.m_z)*mi2;
	cs=CS_Parameters(kt2,ii.m_z,ii.m_y,ii.m_phi,ii.m_z,Q2,3,kin);
      }
    }
  }
  cs.m_col=col;
  KernelWeight(i,j,k,mo,cs);
  return cs;
}

double CS_Cluster_Definitions::GetX
(const Cluster_Leg *l,Splitting_Function_Base *const sf) const
{
  const Vec4D &p(l->Mom());
  if (p.PPlus()<p.PMinus()) {
    if (sf) sf->Lorentz()->SetBeam(0);
    return -p.PPlus()/rpa->gen.PBeam(0).PPlus();
  }
  if (sf) sf->Lorentz()->SetBeam(1);
  return -p.PMinus()/rpa->gen.PBeam(1).PMinus();
}

Flavour CS_Cluster_Definitions::ProperFlav(const Flavour &fl) const
{
  Flavour pfl(fl);
  switch (pfl.Kfcode()) {
  case kf_gluon_qgc: pfl=Flavour(kf_gluon); break;
  case kf_h0_qsc: pfl=Flavour(kf_h0); break;
  case kf_Z_qgc: pfl=Flavour(kf_Z); break;
  case kf_Wplus_qgc: pfl=Flavour(kf_Wplus);
    if (fl.IsAnti()) pfl=pfl.Bar(); break;
  default: break;
  }
  return pfl;
}

void CS_Cluster_Definitions::KernelWeight
(const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,
 CS_Parameters &cs) const
{
  const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFFMap());
  if (cs.m_mode==2) cmap=&p_shower->GetSudakov()->FFIMap();
  else if (cs.m_mode==1) cmap=cs.m_col>0?&p_shower->GetSudakov()->IFFMap():
			   &p_shower->GetSudakov()->FIFMap();
  else if (cs.m_mode==3) cmap=cs.m_col>0?&p_shower->GetSudakov()->IFIMap():
			   &p_shower->GetSudakov()->FIIMap();
  SF_EEE_Map::const_iterator eees(cmap->find(ProperFlav(i->Flav())));
  if (eees==cmap->end()) {
    msg_Debugging()<<"No splitting function, skip kernel weight calc for "
                   <<ProperFlav(i->Flav())<<"("<<i->Flav()<<").\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  SF_EE_Map::const_iterator ees(eees->second.find(ProperFlav(j->Flav())));
  if (ees==eees->second.end()) {
    msg_Debugging()<<"No splitting function, skip kernel weight calc for ["
                   <<ProperFlav(i->Flav())<<"("<<i->Flav()<<"), "
                   <<ProperFlav(j->Flav())<<"("<<j->Flav()<<")].\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  SF_E_Map::const_iterator es(ees->second.find(ProperFlav(mo)));
  if (es==ees->second.end()) {
    msg_Debugging()<<"No splitting function, skip kernel weight calc for ["
                   <<ProperFlav(i->Flav())<<"("<<i->Flav()<<"), "
                   <<ProperFlav(j->Flav())<<"("<<j->Flav()<<")] ["
                   <<ProperFlav(mo)<<"("<<mo<<")].\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  Splitting_Function_Base *cdip(es->second);
  Flavour fls((k->Id()&3)?ProperFlav(k->Flav()).Bar():ProperFlav(k->Flav()));
  if (!cdip->Coupling()->AllowSpec(fls)) {
    msg_Debugging()<<"Invalid spectator "<<fls<<"\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  double Q2=dabs((i->Mom()+j->Mom()+k->Mom()).Abs2());
  cs.p_sf=cdip;
  p_shower->SetMS(p_ms);
  cdip->SetFlavourSpec(fls);
  cs.m_mu2=Max(cs.m_kt2,cs.m_mode&1?
	       p_shower->GetSudakov()->ISPT2Min():
	       p_shower->GetSudakov()->FSPT2Min());
  cs.m_idi=i->Id();
  cs.m_idj=j->Id();
  cs.m_idk=k->Id();
  if (!(m_mode&1)) return;
  double scale=cs.m_kt2, eta=1.0;
  if (cs.m_mode==1) eta=GetX(i,cdip)*cs.m_z;
  else if (cs.m_mode==2) eta=GetX(k,cdip)*(1.0-cs.m_y);
  else if (cs.m_mode==3) eta=GetX(i,cdip)*cs.m_z;
  Color_Info ci(i->Col(),j->Col(),k->Col());
  cs.m_wk=(*cdip)(cs.m_z,cs.m_y,eta,scale,Q2,ci);
  if (cs.m_wk<=0.0 || IsBad(cs.m_wk) || 
      (m_amode==1 && !cdip->On()))
    cs.m_wk=sqrt(std::numeric_limits<double>::min());
  cs.m_ws=cs.m_kt2/cs.m_wk;
  msg_Debugging()<<"Kernel weight [m="
		 <<cs.m_mode<<",c="<<cs.m_col<<"] ( x = "<<eta
		 <<" ) {\n  "<<*i<<"\n  "<<*j<<"\n  "<<*k
		 <<"\n} -> w = "<<cs.m_wk<<" ("<<cs.m_ws<<")\n";
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin,const int kmode)
{
  p_ms=ms;
  if (i>j) std::swap<int>(i,j);
  Vec4D_Vector after(ampl.Legs().size()-1);
  double mi2=p_ms->Mass2(ampl.Leg(i)->Flav());
  double mj2=p_ms->Mass2(ampl.Leg(j)->Flav());
  double mk2=p_ms->Mass2(ampl.Leg(k)->Flav()), mij2=p_ms->Mass2(mo);
  double mb2=i<2?p_ms->Mass2(ampl.Leg(1-i)->Flav()):0.0;
  Vec4D pi(ampl.Leg(i)->Mom()), pj(ampl.Leg(j)->Mom());
  Vec4D pk(ampl.Leg(k)->Mom()), pb(i<2?ampl.Leg(1-i)->Mom():Vec4D());
  if (i>1 && mi2>10.0 && !ampl.Leg(i)->Flav().Strong()) mi2=pi.Abs2();
  if (j>1 && mj2>10.0 && !ampl.Leg(j)->Flav().Strong()) mj2=pj.Abs2();
  if (k>1 && mk2>10.0 && !ampl.Leg(k)->Flav().Strong()) mk2=pk.Abs2();
  bool sk(true);
  if (i>1 && j>1 && mij2>10.0 && !mo.Strong()) {
    mij2=(pi+pj).Abs2();
    pk[0]=pk[0]<0.0?-pk.PSpat():pk.PSpat();
    if (mk2) sk=false;
    mk2=0.0;
  }
  Kin_Args lt;
  if (i>1) {
    if (k>1) lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,2|(kin?4:0));
    else lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,2|(kin?4:0));
    if (k<=1) {
      Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
      if (lt.m_pk.PPlus()>sum.PPlus() ||
	  lt.m_pk.PMinus()>sum.PMinus()) return Vec4D_Vector();
    }
  }
  else {
    if (k>1) lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-pb,2|(kin?4:0));
    else lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,2|(kin?4:0));
    Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
    if (lt.m_pi.PPlus()>sum.PPlus() ||
	lt.m_pi.PMinus()>sum.PMinus()) return Vec4D_Vector();
  }
  if (lt.m_stat<0) return Vec4D_Vector();
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)j) continue;
    if (m==(size_t)i) after[l]=i>1?lt.m_pi:-lt.m_pi;
    else if (m==(size_t)k && sk) after[l]=k>1?lt.m_pk:-lt.m_pk;
    else after[l]=lt.m_lam*ampl.Leg(m)->Mom();
    ++l;
  }
  return after;
}

namespace MCATNLO {

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs)
  {
    return str<<"CS{kt="<<sqrt(cs.m_kt2)<<",z="<<cs.m_z<<",phi="<<cs.m_phi
	      <<",mode="<<cs.m_mode<<",kin="<<cs.m_kin<<"}";
  }

}
