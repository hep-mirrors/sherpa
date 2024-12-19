#include "ALARIC++/Shower/Cluster_Definitions.H"

#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Phys/Recoil_Definition.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"

#include <numeric>

using namespace ALARIC;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

const double s_uxeps=1.0e-3;

Cluster_Definitions::Cluster_Definitions(Shower *const shower):
  p_shower(shower), m_mode(0), m_amode(0) {}

Cluster_Param Cluster_Definitions::Cluster(const Cluster_Config &ca)
{
  DEBUG_FUNC(ca);
  p_ms=ca.p_ms;
  int i(ca.m_i), j(ca.m_j), swap(j<ca.p_ampl->NIn() && j<i);
  if (swap) std::swap<int>(i,j);
  double ws, mu2;
  int type((i<ca.p_ampl->NIn()?1:0)|(ca.m_k<ca.p_ampl->NIn()?2:0));
  Splitting s(KT2(*ca.p_ampl,i,j,ca.m_k,ca.m_mo,ca.m_n,
		  type,(swap?2:0)|(ca.m_mode<<2),ws,mu2));
  bool iss=i<ca.p_ampl->NIn() || j<ca.p_ampl->NIn();
  if (s.m_t>0.0 && s.m_t<sqr(rpa->gen.Ecms()))
    return Cluster_Param(this,ws,s.m_t,mu2,0,ca.m_n,0,s.m_p);
  return Cluster_Param(this,0.0,-1.0);
}

Splitting Cluster_Definitions::KT2
(const ATOOLS::Cluster_Amplitude &ampl,
 int i,int j,int k,const ATOOLS::Flavour &mo,const int id,
 const int type,const int mode,double &ws,double &mu2)
{
  const ATOOLS::Cluster_Leg *li=ampl.Leg(i), *lj=ampl.Leg(j), *lk=ampl.Leg(k);
  Parton c(NULL,li->Flav(),li->Mom(),Color(li->Col().m_i,li->Col().m_j));
  Parton s(NULL,lk->Flav(),lk->Mom(),Color(lk->Col().m_i,lk->Col().m_j));
  Parton n(NULL,lj->Flav(),lj->Mom(),Color(lj->Col().m_i,lj->Col().m_j));
  if (type&1) c.SetBeam(li->Id()&3);
  if (type&2) s.SetBeam(lk->Id()&3);
  Splitting sp(&c,&s);
  sp.SetType();
  sp.m_rcl=p_shower->Recoil()->RecoilTags(&ampl);
  for (size_t l(0);l<ampl.Legs().size();++l)
    if (sp.m_rcl[l]&2) sp.m_Kt+=ampl.Leg(l)->Mom();
  sp.m_eta=c.GetXB();
  sp.p_n=&n;
  sp.m_clu=1;
  sp.m_cpl=p_shower->CouplingScheme();
  sp.m_kfac=(mode&(16<<2))?0:p_shower->KFactorScheme();
  sp.m_pi=li->Mom();
  sp.m_pj=lj->Mom();
  sp.m_pk=lk->Mom();
  Kernel *sk(NULL);
  const Kernel_Vector &sks(p_shower->GetKernels(mo));
  for (size_t n(0), i(0);i<sks.size();++i)
    if ((sks[i]->Type()&1)==(type&1) &&
	sks[i]->LF()->Flav(1)==li->Flav() &&
	sks[i]->LF()->Flav(2)==lj->Flav() &&
	sks[i]->GF()->Allowed(sp)) {
      if (n++==id) {
	sp.p_sk=sk=sks[i];
	break;
      }
    }
  if (sk==NULL) return Splitting(NULL,NULL,-1.0);
  ws=0.0;
  sk->LF()->SetMS(ampl.MS());
  if (!sk->LF()->SetLimits(sp) ||
      !sk->LF()->Cluster(sp,1|2,i,j,k,&ampl)) {
    sp.m_t=-1.0;
    return sp;
  }
  msg_Debugging()<<"c/s mom: "<<sp.m_p[i<j?i:j]
		 <<" "<<sp.m_p[(i<j?j:i)<k?k-1:k]<<"\n";
  msg_Debugging()<<"Mom sum: "<<std::accumulate
    (sp.m_p.begin(),sp.m_p.end(),Vec4D())<<"\n";
  msg_Debugging()<<"Splitting: t = "<<sp.m_t<<" = "<<sqrt(sp.m_t)
		 <<" ^ 2, z = "<<sp.m_z<<", phi = "<<sp.m_phi<<"\n";
  c.SetMom(sp.m_p[i<j?i:j]);
  s.SetMom(sp.m_p[(i<j?j:i)<k?k-1:k]);
  ws=sk->Value(sp);
  mu2=sk->GF()->TrueScale(sp);
  msg_Debugging()<<"Scale: "<<sqrt(mu2)<<" <- "
		 <<sqrt(sk->GF()->Scale(sp))<<"\n";
  msg_Debugging()<<"Kernel: "<<ws<<" ( kfac = "<<sp.m_kfac
		 <<" )  <-  "<<sk->Class()<<"\n";
  if (p_shower->KFactorScheme() &&
      sp.m_t>p_shower->TMin(type&1)) {
    sp.m_kfac=0;
    double K=ws/sk->Value(sp);
    msg_Debugging()<<"     K: "<<K<<" ( kfac = "<<sp.m_kfac<<" )\n";
    if (K>0.0 && !IsEqual(K,1.0)) {
      sp.m_clu=0;
      mu2=sk->GF()->Solve(K*sk->GF()->Coupling(sp));
    }
  }
  if (ws) ws/=sp.m_t;
  return sp;
}

Flavour Cluster_Definitions::ProperFlav(const Flavour &fl) const
{
  Flavour pfl(fl);
  switch (pfl.Kfcode()) {
  case kf_gluon_qgc: pfl=Flavour(kf_gluon); break;
  default: break;
  }
  return pfl;
}
