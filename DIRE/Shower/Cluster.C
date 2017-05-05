#include "DIRE/Shower/Cluster.H"

#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace DIRE;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

const double s_uxeps=1.0e-3;

Cluster::Cluster(Shower *const shower):
  p_shower(shower), m_mode(0), m_amode(0) {}

CParam Cluster::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin,const int kmode)
{
  DEBUG_FUNC("i = "<<i<<", j = "<<j<<", k = "<<k);
  p_shower->SetMS(p_ms=ms);
  int swap(j<ampl.NIn() && j<i);
  if (swap) std::swap<int>(i,j);
  double ws, mu2;
  int type((i<ampl.NIn()?1:0)|(k<ampl.NIn()?2:0));
  Splitting s(KT2(ampl.Leg(i),ampl.Leg(j),ampl.Leg(k),mo,
		  kin,type,swap?2:0,ws,mu2));
  if (s.m_t) return CParam(s.m_t,ws,s.m_eta,mu2,s.m_kin,0);
  double t(dabs((ampl.Leg(i)->Mom()+ampl.Leg(j)->Mom()).Abs2()));
  return CParam(t,t,0.0,t,0,0);
}

Splitting Cluster::KT2
(const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
 const ATOOLS::Cluster_Leg *lk,const ATOOLS::Flavour &mo,
 const int kin,const int type,const int mode,double &ws,double &mu2)
{
  Parton c(NULL,li->Flav(),li->Mom(),Color(li->Col().m_i,li->Col().m_j));
  Parton s(NULL,lk->Flav(),lk->Mom(),Color(lk->Col().m_i,lk->Col().m_j));
  Parton n(NULL,lj->Flav(),lj->Mom(),Color(lj->Col().m_i,lj->Col().m_j));
  if (type&1) c.SetBeam(li->Id()&3);
  if (type&2) s.SetBeam(lk->Id()&3);
  Splitting sp(&c,&s);
  sp.m_eta=c.GetXB();
  sp.p_n=&n;
  sp.m_kin=kin;
  sp.m_clu=1;
  sp.m_type=type;
  sp.m_cpl=p_shower->CouplingScheme();
  Kernel *sk(p_shower->GetKernel(sp,(mode&2)?1:0));
  if (sk==NULL) return Splitting();
  sk->LF()->Cluster(sp,mode&~2);
  msg_Debugging()<<"Splitting: t = "<<sp.m_t<<" = "<<sqrt(sp.m_t)
		 <<" ^ 2, z = "<<sp.m_z<<", phi = "<<sp.m_phi<<"\n"; 
  ws=sk->Value(sp);
  msg_Debugging()<<"Kernel: "<<ws<<" ( kfac = "<<sp.m_kfac<<" )\n";
  sp.m_clu=0;
  mu2=sk->GF()->Scale(sp);
  if (p_shower->KFactorScheme() &&
      mu2>p_shower->TMin(sp.m_type&1)) {
    sp.m_kfac=p_shower->KFactorScheme();
    double K=sk->GF()->K(sp);
    msg_Debugging()<<"     K: "<<K<<" ( kfac = "<<sp.m_kfac<<" )\n";
    sp.m_kfac=0;
    if (!IsZero(K)) mu2=sk->GF()->Solve((1.0+K)*sk->GF()->Coupling(sp));
  }
  if (ws) ws=sp.m_t/ws; 
  else ws=std::numeric_limits<double>::max();
  if (ws<0.0) ws=-ws;
  return sp;
}

Flavour Cluster::ProperFlav(const Flavour &fl) const
{
  Flavour pfl(fl);
  switch (pfl.Kfcode()) {
  case kf_gluon_qgc: pfl=Flavour(kf_gluon); break;
  default: break;
  }
  return pfl;
}

ATOOLS::Vec4D_Vector  Cluster::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin,const int kmode)
{
  p_ms=ms;
  if (j<ampl.NIn() && j<i) std::swap<int>(i,j);
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
