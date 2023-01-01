#include "PHASIC++/Scales/Cluster_Definitions.H"

#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Cluster_Param Cluster_Definitions::Cluster(const Cluster_Config &cc)
{
  Cluster_Param cp(this,0.0);
  const Cluster_Leg *lk(cc.p_ampl->Leg(cc.m_k));
  const Cluster_Leg *li(cc.p_ampl->Leg(cc.m_i)), *lj(cc.p_ampl->Leg(cc.m_j));
  int cpl((li->Flav().Strong()&&lj->Flav().Strong()&&cc.m_mo.Strong())?0:1);
  if ((li->Flav().IsGluon()&&lj->Flav().IsGluon()&&!cc.m_mo.Strong()) ||
      (li->Flav().IsGluon()&&cc.m_mo.IsGluon()&&!lj->Flav().Strong()) ||
      (lj->Flav().IsGluon()&&cc.m_mo.IsGluon()&&!li->Flav().Strong())) cpl=2;
  if (cc.m_j<cc.m_i) std::swap<const Cluster_Leg*>(li,lj);
  Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom());
  double mi2=cc.p_ms->Mass2(li->Flav()), mj2=cc.p_ms->Mass2(lj->Flav());
  double mk2=cc.p_ms->Mass2(lk->Flav()), mij2=cc.p_ms->Mass2(cc.m_mo);
  double z=-1.0, t=2.0*dabs(pi*pj), mu=t, Q2=2.0*dabs(pi*pj+pk*(pi+pj));
  if (li->Stat()&4) mi2=pi.Abs2();
  if (lj->Stat()&4) mj2=pj.Abs2();
  if (lk->Stat()&4) mk2=pk.Abs2();
  int ib(cc.m_i<cc.m_j?cc.m_i:cc.m_j), kin(1);
  if (cpl==1 && cc.m_mo.Width()) {
    mij2=(pi+pj).Abs2();
    kin=0;
  }
  if ((li->Id()&3)==0) {
    if ((lj->Id()&3)==0) {
      t+=mi2+mj2;
      if ((lk->Id()&3)==0) {
	Kin_Args ffp(ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,3));
	if (cpl==0) {
	  if (m_pstype) {
	    t=Q2*ffp.m_y*(1.0-ffp.m_y)*(cc.m_i<cc.m_j?(1.0-ffp.m_z):ffp.m_z);
	    z=1.0-(cc.m_i<cc.m_j?(1.0-ffp.m_z):ffp.m_z)*(1.0-ffp.m_y);
	  }
	  else {
	    t=Q2*ffp.m_y*ffp.m_z*(1.0-ffp.m_z);
	    z=(cc.m_i<cc.m_j?ffp.m_z:1.0-ffp.m_z);
	  }
	}
	if (ffp.m_stat<0) t=-1.0;
	cp=Cluster_Param(this,0.0,t,t,cpl,kin,0,ffp.m_pi,ffp.m_pk,ffp.m_lam);
      }
      else {
	Kin_Args fip(ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,3));
	if (cpl==0) {
	  if (m_pstype) {
	    t=Q2*fip.m_y*(cc.m_i<cc.m_j?1.0-fip.m_z:fip.m_z);
	    z=(cc.m_i<cc.m_j?fip.m_z:1.0-fip.m_z);
	  }
	  else {
	    t=Q2*fip.m_y/(1.0-fip.m_y)*fip.m_z*(1.0-fip.m_z);
	    z=(cc.m_i<cc.m_j?fip.m_z:1.0-fip.m_z);
	  }
	}
	if (fip.m_pk[0]<0.0 || fip.m_stat<0 ||
	    fip.m_pk[0]>rpa->gen.PBunch(cc.m_k)[0]) t=-1.0;
	cp=Cluster_Param(this,0.0,t,t,cpl,kin,0,fip.m_pi,-fip.m_pk,fip.m_lam);
      }
    }
  }
  else {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	const Vec4D &pb(cc.p_ampl->Leg(cc.m_j<cc.m_i?1-cc.m_j:1-cc.m_i)->Mom());
	Kin_Args ifp(ClusterIFDipole(mi2,mj2,mij2,mk2,0.0,-pi,pj,pk,-pb,3|(kin?4:0)));
	if (cpl==0) {
	  if (m_pstype) {
	    mu=(t=Q2*ifp.m_y*(1.0-ifp.m_z))/ifp.m_z;
	    z=ifp.m_z;
	  }
	  else {
	    mu=t=Q2*ifp.m_y/ifp.m_z*(1.0-ifp.m_z);
	    z=ifp.m_z;
	  }
	}
	if (ifp.m_pi[0]<0.0 || ifp.m_stat<0 ||
	    ifp.m_pi[0]>rpa->gen.PBunch(ib)[0]) t=-1.0;
	cp=Cluster_Param(this,0.0,t,mu,cpl,kin,0,-ifp.m_pi,ifp.m_pk,ifp.m_lam);
      }
      else {
	Kin_Args iip(ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,3|(kin?4:0)));
	if (cpl==0) {
	  if (m_pstype) {
	    mu=(t=Q2*iip.m_y*(1.0-iip.m_z-iip.m_y))/iip.m_z;
	    z=iip.m_z+iip.m_y;
	  }
	  else {
	    mu=t=Q2*iip.m_y/iip.m_z*(1.0-iip.m_z);
	    z=iip.m_z;
	  }
	}
	if (iip.m_pi[0]<0.0 || iip.m_stat<0 ||
	    iip.m_pi[0]>rpa->gen.PBunch(ib)[0]) t=-1.0;
	cp=Cluster_Param(this,0.0,t,mu,cpl,kin,0,-iip.m_pi,-iip.m_pk,iip.m_lam);
      }
    }
  }
  if (cp.m_kt2<=0.0) return cp;
  cp.m_op=(8.0*M_PI)/dabs((pi+pj).Abs2()-cc.p_ms->Mass2(cc.m_mo));
  int isi(cc.m_i<cc.p_ampl->NIn()||cc.m_j<cc.p_ampl->NIn());
  if (cp.m_cpl==1) {
    cp.m_op*=(*MODEL::aqed)(t);
    cp.m_stat=4;
  }
  else {
    cp.m_op*=(*MODEL::as)(sqr(rpa->gen.Ecms()));
    if (m_kfac && !m_nproc &&
	((isi || cc.m_i<cc.m_j)?lj:li)->Flav().IsGluon() &&
	!(isi && cc.m_j<cc.m_i &&
	  li->Flav().IsGluon() && lj->Flav().IsGluon())) {
      double nf=MODEL::as->Nf(cp.m_kt2), B0=11.0/2.0-1.0/3.0*nf;
      double G2=(67.0/6.0-sqr(M_PI)/2.0)-5.0/9.0*nf;
      cp.m_op*=1.0+(*MODEL::as)(cp.m_kt2)/(2.0*M_PI)*G2;
      cp.m_mu2*=exp(-G2/B0);
    }
    if (abs(cc.m_mo.StrongCharge())==3) {
      if (!isi) cp.m_op*=4.0/3.0;
      else cp.m_op*=li->Flav().StrongCharge()==8?0.5:4.0/3.0;
    }
    else {
      if (li->Flav().StrongCharge()==8) cp.m_op*=3.0;
      else cp.m_op*=0.5*(isi?4.0/3.0:0.5);
    }
  }
  if (!isi && cc.m_j<cc.m_i) std::swap<const Cluster_Leg*>(li,lj);
  if (cc.m_mo.IsVector()) {
    if (li->Flav().IsVector()) {
      if (lj->Flav().IsVector()) {// VVV
	if (!isi) cp.m_op*=(1.0-z)/(sqr(1.0-z)+cp.m_kt2/Q2)-1.0+z*(1.0-z)/2.0;
	else {
	  if (cc.m_i<cc.m_j) cp.m_op*=(1.0-z)/(sqr(1.0-z)+cp.m_kt2/Q2)-1.0
			       +(z/(z*z+cp.m_kt2/Q2)-1.0)/2.0;
	  else cp.m_op*=(z/(z*z+cp.m_kt2/Q2)-1.0)/2.0+z*(1.0-z);
	}
      }
    }
    else if (li->Flav().IsFermion()) {
      if (lj->Flav().IsFermion()) {
	if (!isi) cp.m_op*=1.0-2.0*z*(1.0-z);// VFF
	else {// FVF
	  if (isi && ((cc.m_i<cc.m_j)^li->Flav().IsAnti())) cp.m_op=0.0;
	  else cp.m_op*=2.0*z/(z*z+cp.m_kt2/Q2)-2.0+z;
	}
      }
    }
  }
  else if (cc.m_mo.IsFermion()) {
    if (li->Flav().IsVector()) {
      if (lj->Flav().IsFermion()) {
	if (!isi) cp.m_op=0.0;// FVF
	else {// VFF
	  if (isi && ((cc.m_i<cc.m_j)^lj->Flav().IsAnti())) cp.m_op=0.0;
	  else cp.m_op*=1.0-2.0*z*(1.0-z);
	}
      }
    }
    else if (li->Flav().IsFermion()) {
      if (lj->Flav().IsVector()) {// FFV
	if (isi && ((cc.m_i>cc.m_j)^li->Flav().IsAnti())) cp.m_op=0.0;
	else cp.m_op*=2.0*(1.0-z)/(sqr(1-z)+cp.m_kt2/Q2)-2.0+(1.0-z);
      }
    }
  }
  return cp;
}
