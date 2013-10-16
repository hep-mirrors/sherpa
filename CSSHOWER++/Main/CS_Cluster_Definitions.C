#include "CSSHOWER++/Main/CS_Cluster_Definitions.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

CS_Cluster_Definitions::CS_Cluster_Definitions
(Shower *const shower,const int kmode,const int meweight):
  p_shower(shower), m_kmode(kmode), m_meweight(meweight) {}

CParam CS_Cluster_Definitions::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin,const int mode)
{
  DEBUG_FUNC("kin = "<<kin<<", mode = "<<mode);
  CS_Parameters cs(KT2(&ampl,ampl.Leg(i),ampl.Leg(j),ampl.Leg(k),mo,ms,kin,mode|m_kmode));
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
  double Q2=(pi+pj+pk).Abs2();
  double mb2=p_b->Mom().Abs2(), mfb2=p_ms->Mass2(p_b->Flav());
  if ((mfb2==0.0 && IsZero(mb2,1.0e-6)) || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
  double mi2=pi.Abs2(), mfi2=p_ms->Mass2(i->Flav());
  double mj2=pj.Abs2(), mfj2=p_ms->Mass2(j->Flav());
  double mk2=pk.Abs2(), mfk2=p_ms->Mass2(k->Flav());
  if ((mfi2==0.0 && IsZero(mi2,1.0e-6)) || IsEqual(mi2,mfi2,1.0e-6)) mi2=mfi2;
  if ((mfj2==0.0 && IsZero(mj2,1.0e-6)) || IsEqual(mj2,mfj2,1.0e-6)) mj2=mfj2;
  if ((mfk2==0.0 && IsZero(mk2,1.0e-6)) || IsEqual(mk2,mfk2,1.0e-6)) mk2=mfk2;
  double mij2=p_ms->Mass2(mo);
  if (kmode&1) {
    mij2=(pi+pj).Abs2();
    pk*=(pi+pj).PSpat()/pk.PSpat();
    pk[0]=pk[0]<0.0?-pk.PSpat():pk.PSpat();
    if (IsZero(Q2,1.0e-3)) pk=pk*4.0;
    mk2=0.0;
    kin=0;
  }
  Q2=(pi+pj+pk).Abs2();
  Kin_Args lt;
  CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),
		   1.0,1.0,0.0,0.0,0.0,
		   ((i->Id()&3)?1:0)|((k->Id()&3)?2:0),kin,kmode&1);
  cs.m_wk=cs.m_ws=-1.0;
  if ((i->Id()&3)==0) {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) {
	lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,1|(kin?4:0));
	if (lt.m_stat!=1) return cs;
	double kt2=p_shower->KinFF()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,1.0,Q2,0,kin,kmode&1);
	cs.m_pk=pk;
      }
      else {
	lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,1|8|(kin?4:0));
	Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
	if ((k==ampl->Leg(0) && lt.m_pk[3]<0.0) ||
	    (k==ampl->Leg(1) && lt.m_pk[3]>0.0) ||
	    lt.m_pk[0]<0.0 || lt.m_y>1.0 || lt.m_stat!=1) return cs;
	double kt2=p_shower->KinFI()->GetKT2(Q2,1.0-lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,1.0-lt.m_y,Q2,2,kin,kmode&1);
	cs.m_pk=-pk;
      }
    }
  }
  else {
    if ((j->Id()&3)==0) {
      Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
      if ((k->Id()&3)==0) {
	lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-p_b->Mom(),3|(kin?4:0));
	if ((kmode&1) && lt.m_mode) lt.m_stat=-1;
	if ((i==ampl->Leg(0) && lt.m_pi[3]<0.0) ||
	    (i==ampl->Leg(1) && lt.m_pi[3]>0.0) ||
	    lt.m_pi[0]<0.0 || lt.m_z<0.0 || lt.m_stat!=1) return cs;
	double kt2=p_shower->KinIF()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,lt.m_z,Q2,1,lt.m_mode,kmode&1);
	cs.m_pk=pk;
      }
      else {
	lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,3|(kin?4:0));
	if ((i==ampl->Leg(0) && lt.m_pi[3]<0.0) ||
	    (i==ampl->Leg(1) && lt.m_pi[3]>0.0) ||
	    lt.m_pi[0]<0.0 || lt.m_z<0.0 || lt.m_stat!=1) return cs;
	double kt2=p_shower->KinII()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,lt.m_z,Q2,3,kin,kmode&1);
	cs.m_pk=-pk;
      }
    }
  }
  cs.m_col=col;
  KernelWeight(i,j,k,mo,cs,kmode);
  if (cs.m_wk>0.0 &&
      (((m_meweight&1) && (kmode&2)) || (kmode&4))) {
    Cluster_Amplitude *campl(Cluster_Amplitude::New());
    campl->SetProcs(ampl->Procs<void>());
    campl->Decays()=ampl->Decays();
    campl->SetNIn(ampl->NIn());
    campl->SetMuR2(rpa->gen.CplScale());
    campl->SetMuF2(rpa->gen.CplScale());
    campl->SetQ2(rpa->gen.CplScale());
    for (size_t l(0), m(0);m<ampl->Legs().size();++m) {
      Cluster_Leg *lm(ampl->Leg(m));
      if (lm==j) continue;
      if (lm==i) campl->CreateLeg((i->Id()&3)?-lt.m_pi:lt.m_pi,mo,CombineColors(i,j,k,mo));
      else if (lm==k) campl->CreateLeg((k->Id()&3)?-lt.m_pk:lt.m_pk,lm->Flav(),lm->Col());
      else campl->CreateLeg(lt.m_lam*ampl->Leg(m)->Mom(),lm->Flav(),lm->Col());
      ++l;
    }
#ifndef DEBUG__Differential
    int olv(msg->Level());
    msg->SetLevel(2);
#endif
    double me=Differential(campl,kmode);
#ifndef DEBUG__Differential
    msg->SetLevel(olv);
#endif
    if (me) cs.m_ws/=me;
    else cs.m_ws=-1.0;
    campl->Delete();
    msg_Debugging()<<"ME = "<<me<<" -> "<<1.0/cs.m_ws<<" ("
		   <<(cs.m_ws>0.0?sqrt(cs.m_ws):-sqrt(-cs.m_ws))<<")\n";
  }
  return cs;
}

double CS_Cluster_Definitions::Differential
(Cluster_Amplitude *const ampl,const int kmode) const
{
  NLOTypeStringProcessMap_Map *procs
    (ampl->Procs<NLOTypeStringProcessMap_Map>());
  if (procs==NULL) return 1.0;
  nlo_type::code type=nlo_type::lo;
  if (procs->find(type)==procs->end()) return 0.0;
  Process_Base::SortFlavours(ampl);
  int rm(ampl->Leg(0)->Mom()[3]<0.0?0:1024);
  std::string pname(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit((*(*procs)[type]).find(pname));
  if (pit==(*(*procs)[type]).end()) return 0.0;
  int cm(pit->second->NOut()>pit->second->Info().m_fi.NMinExternal()?1|64:0);
  if (kmode&4) cm&=~1;
  double meps=pit->second->Differential(*ampl,cm|2|4|((m_meweight&2)?64:0)|rm);
  meps*=pit->second->SymFac();
  return meps;
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
 CS_Parameters &cs,const int kmode) const
{
  const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFFMap());
  if (cs.m_mode==2) cmap=&p_shower->GetSudakov()->FFIMap();
  else if (cs.m_mode==1) cmap=cs.m_col>0?&p_shower->GetSudakov()->IFFMap():
			   &p_shower->GetSudakov()->FIFMap();
  else if (cs.m_mode==3) cmap=cs.m_col>0?&p_shower->GetSudakov()->IFIMap():
			   &p_shower->GetSudakov()->FIIMap();
  SF_EEE_Map::const_iterator eees(cmap->find(ProperFlav(i->Flav())));
  if (eees==cmap->end()) {
    msg_Debugging()<<"No splitting function (i), skip kernel weight calc\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  SF_EE_Map::const_iterator ees(eees->second.find(ProperFlav(j->Flav())));
  if (ees==eees->second.end()) {
    msg_Debugging()<<"No splitting function (j), skip kernel weight calc\n";
    cs.m_ws=cs.m_wk=-1.0;
    return;
  }
  SF_E_Map::const_iterator es(ees->second.find(ProperFlav(mo)));
  if (es==ees->second.end()) {
    msg_Debugging()<<"No splitting function (k), skip kernel weight calc\n";
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
  cs.m_mu2*=cdip->Coupling()->CplFac(cs.m_mu2);
  if (!cdip->On()) cs.m_mu2=Max(cs.m_mu2,sqr(mo.Mass()));
  if (!(kmode&2)) return;
  if (!cdip->On()) {
    if (AMode()==1) cs.m_wk=-1.0;
    else cs.m_wk=sqrt(std::numeric_limits<double>::min());
    cs.m_ws=1.0/cs.m_wk;
    msg_Debugging()<<"No Kernel. Set weight "<<cs.m_ws<<".\n";
  }
  else {
  double scale=cs.m_kt2, eta=1.0;
  if (cs.m_mode==1) eta=GetX(i,cdip)*cs.m_z;
  else if (cs.m_mode==2) eta=GetX(k,cdip)*(1.0-cs.m_y);
  else if (cs.m_mode==3) eta=GetX(i,cdip)*cs.m_z;
  cs.m_wk=(*cdip)(cs.m_z,cs.m_y,eta,-1.0,Q2)*
    cdip->MEPSWeight(cs.m_z,cs.m_y,eta,-1.0,Q2);
  cs.m_wk*=cdip->SymFac();
  if (cs.m_wk<=0.0 || IsBad(cs.m_wk))
    cs.m_wk=sqrt(std::numeric_limits<double>::min());
  cs.m_ws=1.0/cs.m_wk;
  msg_Debugging()<<"Kernel weight (A="<<AMode()
		 <<") [m="<<cs.m_mode<<",c="<<cs.m_col<<"] ( x = "<<eta
		 <<" ) {\n  "<<*i<<"\n  "<<*j<<"\n  "<<*k
		 <<"\n} -> w = "<<cs.m_wk<<" ("<<cs.m_ws<<")\n";
  }
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int ikin,const int kmode)
{
  p_ms=ms;
  int kin(ikin);
  if (i>j) std::swap<int>(i,j);
  Vec4D_Vector after(ampl.Legs().size()-1);
  double mb2(0.0);
  if (i<2) {
    mb2=ampl.Leg(1-i)->Mom().Abs2();
    double mfb2(p_ms->Mass2(ampl.Leg(1-i)->Flav()));
    if ((mfb2==0.0 && IsZero(mb2,1.0e-6)) || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
  }
  Vec4D pi(ampl.Leg(i)->Mom()), pj(ampl.Leg(j)->Mom());
  Vec4D pk(ampl.Leg(k)->Mom()), pb(i<2?ampl.Leg(1-i)->Mom():Vec4D());
  double mi2=pi.Abs2(), mfi2=p_ms->Mass2(ampl.Leg(i)->Flav());
  double mj2=pj.Abs2(), mfj2=p_ms->Mass2(ampl.Leg(j)->Flav());
  double mk2=pk.Abs2(), mfk2=p_ms->Mass2(ampl.Leg(k)->Flav());
  if ((mfi2==0.0 && IsZero(mi2,1.0e-6)) || IsEqual(mi2,mfi2,1.0e-6)) mi2=mfi2;
  if ((mfj2==0.0 && IsZero(mj2,1.0e-6)) || IsEqual(mj2,mfj2,1.0e-6)) mj2=mfj2;
  if ((mfk2==0.0 && IsZero(mk2,1.0e-6)) || IsEqual(mk2,mfk2,1.0e-6)) mk2=mfk2;
  double mij2=p_ms->Mass2(mo);
  bool sk(true);
  if (kmode&1) {
    mij2=(pi+pj).Abs2();
    pk[0]=pk[0]<0.0?-pk.PSpat():pk.PSpat();
    if (mk2) sk=false;
    mk2=0.0;
    kin=0;
  }
  Kin_Args lt;
  if (i>1) {
    if (k>1) lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,2|(kin?4:0));
    else lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,2|(kin?4:0));
    if ((k==0 && lt.m_pk[3]<0.0) ||
	(k==1 && lt.m_pk[3]>0.0) || lt.m_pk[0]<0.0) return Vec4D_Vector();
  }
  else {
    if (k>1) {
      lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-pb,2|(kin?4:0));
      if ((kmode&1) && lt.m_mode) lt.m_stat=-1;
    }
    else lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,2|(kin?4:0));
    if ((i==0 && lt.m_pi[3]<0.0) ||
	(i==1 && lt.m_pi[3]>0.0) || lt.m_pi[0]<0.0) return Vec4D_Vector();
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

bool CS_Cluster_Definitions::CheckColors
(const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
 const ATOOLS::Cluster_Leg *lk,const ATOOLS::Flavour &mo) const
{
  if (mo.StrongCharge()==8) {
    if (!lk->Flav().Strong()) return false;
  }
  else if (mo.Strong()) {
    if (!(lk->Flav().StrongCharge()==8 ||
	  lk->Flav().StrongCharge()==-mo.StrongCharge())) return false;
  }
  else {
    if (lk->Flav().StrongCharge()==8) return false;
    if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
	lk->Col().m_i==-1) return true;
    ColorID ci(li->Col()), cj(lj->Col());
    if (ci.m_i==cj.m_j && ci.m_j==0 && cj.m_i==0) return true;
    if (ci.m_j==cj.m_i && ci.m_i==0 && cj.m_j==0) return true;
    return false;
  }
  if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
      lk->Col().m_i==-1) return true;
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i<0 && cj.m_i<0 && ck.m_i<0) return true;
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	  (ci.m_i==cj.m_j && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_i==cj.m_j && 
	  (cj.m_i==ck.m_j || ck.Singlet())) return true;
      if ((ci.m_i==ck.m_j || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	  (ci.m_j==cj.m_i && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_j==cj.m_i && 
	  (cj.m_j==ck.m_i || ck.Singlet())) return true;
      if ((ci.m_j==ck.m_i || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lk->Flav().StrongCharge()==0) return false;
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j && 
	  (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	   (ci.m_j==cj.m_i && lk->Flav().StrongCharge()!=8))) 
	return true;
      if (ci.m_j==cj.m_i && 
	  (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	   (ci.m_i==cj.m_j && lk->Flav().StrongCharge()!=8)))
	return true;
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.m_j==cj.m_i &&
	  (ci.m_i==ck.m_j || ck.Singlet())) return true;
      if ((cj.m_i==ck.m_j || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.m_i==cj.m_j &&
	  (ci.m_j==ck.m_i || ck.Singlet())) return true;
      if ((cj.m_j==ck.m_i || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else {
      return false;
    }
  }
  else {
    if (lj->Flav().StrongCharge()==8 ||
	lk->Flav().StrongCharge()==8) {
      return false;
    }
    return true;
  }
  return false;
}

ColorID CS_Cluster_Definitions::CombineColors
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,const ATOOLS::Flavour &mo) const
{
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i==-1 && cj.m_i==-1 && ck.m_i==-1) return ColorID();
  if (!mo.Strong()) return ColorID(0,0);
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      return ColorID(ci.m_i,cj.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(ci.m_i,0);
      return ColorID(cj.m_i,0);
    }
    else {
      return ColorID(ci.m_i,0);
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,ci.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(0,ci.m_j);
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,ci.m_j);
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j) return ColorID(cj.m_i,ci.m_j);
      if (ci.m_j==cj.m_i) return ColorID(ci.m_i,cj.m_j);
      THROW(fatal_error,"Invalid clustering");
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.Singlet()) return ColorID(cj.m_i,0);
      return ColorID(ci.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.Singlet()) return ColorID(0,cj.m_j);
      return ColorID(0,ci.m_j);
    }
    else {
      THROW(fatal_error,"Invalid combination");
    }
  }
  else {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,0);
    }
  }
  return ColorID();
}

namespace CSSHOWER {

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs)
  {
    return str<<"CS{kt="<<sqrt(cs.m_kt2)<<",z="<<cs.m_z<<",phi="<<cs.m_phi
	      <<",mode="<<cs.m_mode<<",kin="<<cs.m_kin
	      <<",kmode="<<cs.m_kmode<<"}";
  }

}
