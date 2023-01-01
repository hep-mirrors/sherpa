#include "CSSHOWER++/Main/CS_Cluster_Definitions.H"

#include "CSSHOWER++/Showers/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

CS_Cluster_Definitions::CS_Cluster_Definitions
(Shower *const shower,const int kmode,const int pdfcheck,const int kfmode):
  p_shower(shower), m_kmode(kmode), m_pdfcheck(pdfcheck), m_kfmode(kfmode) {}

Cluster_Param CS_Cluster_Definitions::Cluster(const Cluster_Config &ca)
{
  DEBUG_FUNC(ca);
  CS_Parameters cs(KT2(ca.p_ampl,
		       ca.p_ampl->Leg(ca.m_i),
		       ca.p_ampl->Leg(ca.m_j),
		       ca.p_ampl->Leg(ca.m_k),
		       ca.m_mo,ca.p_ms,ca.m_kin,ca.m_mode|m_kmode));
  return Cluster_Param(this,cs.m_wk,cs.m_kt2,cs.m_mu2,cs.m_cpl,
		       cs.m_kin,cs.m_kmode,cs.m_pijt,cs.m_pkt,cs.m_lt);
}

Splitting_Function_Base *CS_Cluster_Definitions::GetSF
(const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,CS_Parameters &cs) const
{
  const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFFMap());
  if (cs.m_mode==2) cmap=&p_shower->GetSudakov()->FFIMap();
  else if (cs.m_mode==1) cmap=cs.m_col>=0?&p_shower->GetSudakov()->IFFMap():
			   &p_shower->GetSudakov()->FIFMap();
  else if (cs.m_mode==3) cmap=cs.m_col>=0?&p_shower->GetSudakov()->IFIMap():
			   &p_shower->GetSudakov()->FIIMap();
  SF_EEE_Map::const_iterator eees(cmap->find(ProperFlav(i->Flav(),i->Id()&3)));
  if (eees==cmap->end()) {
    msg_Debugging()<<"No splitting function (i)\n";
    return NULL;
  }
  SF_EE_Map::const_iterator ees(eees->second.find(ProperFlav(j->Flav(),j->Id()&3)));
  if (ees==eees->second.end()) {
    msg_Debugging()<<"No splitting function (j)\n";
    return NULL;
  }
  SF_E_Map::const_iterator es(ees->second.find(ProperFlav(mo,(i->Id()&3)|(j->Id()&3))));
  if (es==ees->second.end()) {
    msg_Debugging()<<"No splitting function (ij)\n";
    return NULL;
  }
  return es->second;
}

CS_Parameters CS_Cluster_Definitions::KT2
(const ATOOLS::Cluster_Amplitude *ampl,
 const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,
 const ATOOLS::Mass_Selector *const ms,const int ikin,
 const int kmode,const int force)
{
  p_ms=ms;
  int kin(ikin<0?p_shower->KinScheme():ikin), col(1);
  if ((i->Id()&3)<(j->Id()&3)) {
    std::swap<const Cluster_Leg*>(i,j);
    col=-1;
  }
  p_b=ampl->Leg(i==ampl->Leg(0)?1:0);
  const Vec4D pi(i->Mom()), pj(j->Mom()), pk(k->Mom());
  const double Q2=(pi+pj+pk).Abs2();
  double mb2=p_b->Mom().Abs2(), mfb2=p_ms->Mass2(p_b->Flav());
  if (mfb2==0.0 || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
  double mi2=pi.Abs2(), mfi2=p_ms->Mass2(i->Flav());
  double mj2=pj.Abs2(), mfj2=p_ms->Mass2(j->Flav());
  double mk2=pk.Abs2(), mfk2=p_ms->Mass2(k->Flav());
  if ((mfi2==0.0 && IsZero(mi2,1.0e-6)) || IsEqual(mi2,mfi2,1.0e-6)) mi2=mfi2;
  if ((mfj2==0.0 && IsZero(mj2,1.0e-6)) || IsEqual(mj2,mfj2,1.0e-6)) mj2=mfj2;
  if ((mfk2==0.0 && IsZero(mk2,1.0e-6)) || IsEqual(mk2,mfk2,1.0e-6)) mk2=mfk2;
  double mij2=p_ms->Mass2(mo);
  if (kmode&1) {
    mij2=(pi+pj).Abs2();
    kin=0;
  }
  Kin_Args lt;
  CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),
		   1.0,1.0,0.0,0.0,0.0,
		   ((i->Id()&3)?1:0)|((k->Id()&3)?2:0),kin,kmode&1);
  cs.m_kt2=-1.0;
  cs.m_wk=0.0;
  cs.m_col=col;
  if ((i->Id()&3)==0) {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) {
	lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,1|(kin?4:0));
	if (lt.m_stat!=1) if (!force) return cs;
	double kt2=p_shower->KinFF()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,1.0,Q2,0,kin,kmode&1);
	cs.m_pkt=cs.m_pk=lt.m_pk;
	cs.m_pijt=lt.m_pi;
      }
      else {
	lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,1|8|(kin?4:0));
	if (lt.m_pk[3]*pk[3]>0.0 ||
	    lt.m_pk[0]<0.0 || lt.m_stat!=1) if (!force) return cs;
	double kt2=p_shower->KinFI()->GetKT2(Q2,1.0-lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,1.0-lt.m_y,Q2,2,kin,kmode&1);
	cs.m_pkt=-(cs.m_pk=lt.m_pk);
	cs.m_pijt=lt.m_pi;
      }
    }
  }
  else {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) {
	lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-p_b->Mom(),3|(kin?4:0));
	if ((kmode&1) && lt.m_mode) lt.m_stat=-1;
	if (lt.m_pi[3]*pi[3]>0.0 ||
	    lt.m_pi[0]<0.0 || lt.m_z<0.0 || lt.m_stat!=1) if (!force) return cs;
	double kt2=p_shower->KinIF()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,lt.m_z,Q2,1,lt.m_mode,kmode&1);
	cs.m_pkt=cs.m_pk=lt.m_pk;
	cs.m_pijt=-lt.m_pi;
	cs.m_lt=lt.m_lam;
      }
      else {
	lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,3|(kin?4:0));
	if (lt.m_pi[3]*pi[3]>0.0 ||
	    lt.m_pi[0]<0.0 || lt.m_z<0.0 || lt.m_stat!=1) if (!force) return cs;
	double kt2=p_shower->KinII()->GetKT2(Q2,lt.m_y,lt.m_z,mi2,mj2,mk2,mo,j->Flav());
	cs=CS_Parameters(kt2,lt.m_z,lt.m_y,lt.m_phi,lt.m_z,Q2,3,kin,kmode&1);
	cs.m_pkt=-(cs.m_pk=lt.m_pk);
	cs.m_pijt=-lt.m_pi;
	cs.m_lt=lt.m_lam;
      }
    }
  }
  cs.m_col=col;
  KernelWeight(i,j,k,mo,cs,kmode);
  return cs;
}

double CS_Cluster_Definitions::GetX(const Cluster_Leg* l,
                                    Splitting_Function_Base* const sf) const
{
  if ((l->Id() & 3) == 0)
    THROW(fatal_error, "Invalid call for parton ID="+ToString(l->Id()));

  int beam((l->Id() & 1) ? 0 : 1);
  if (sf)
    sf->Lorentz()->SetBeam(beam);

  Vec4D mom( (l->Mom()[0] < 0.0) ? -l->Mom() : l->Mom() );
  return p_shower->ISR()->CalcX(mom);
}

Flavour CS_Cluster_Definitions::ProperFlav(const Flavour &fl,const int anti) const
{
  Flavour pfl(fl);
  switch (pfl.Kfcode()) {
  case kf_gluon_qgc: pfl=Flavour(kf_gluon); break;
  default: break;
  }
  return anti?pfl.Bar():pfl;
}

void CS_Cluster_Definitions::KernelWeight
(const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,
 CS_Parameters &cs,const int kmode) const
{
  Splitting_Function_Base *cdip(GetSF(i,j,k,mo,cs));
  if (cdip==NULL || !cdip->On()) {
    cs.m_wk=0.0;
    return;
  }
  if (m_pdfcheck && (cs.m_mode&1)) {
    int beam=i->Id()&1?0:1;
    if (p_shower->ISR()->PDF(beam) &&
	!p_shower->ISR()->PDF(beam)->Contains(mo)) {
      msg_Debugging()<<"Not in PDF: "<<mo<<".\n";
      cs.m_wk=0.0;
      cs.m_kmode=-1;
      return;
    }
  }
  Flavour fls(ProperFlav(k->Flav(),k->Id()&3));
  if (!(kmode&32) && !cdip->Coupling()->AllowSpec(fls,1)) {
    msg_Debugging()<<"Invalid spectator "<<fls<<"\n";
    cs.m_wk=0.0;
    return;
  }
  double Q2=dabs((i->Mom()+j->Mom()+k->Mom()).Abs2());
  if (Q2<(cs.m_mode&1?
	  p_shower->GetSudakov()->ISPT2Min():
	  p_shower->GetSudakov()->FSPT2Min())) {
    msg_Debugging()<<"Small Q2 "<<Q2<<"\n";
    cs.m_wk=0.0;
    return;
  }
  cs.p_sf=cdip;
  p_shower->SetMS(p_ms);
  cdip->SetFlavourSpec(fls);
  cs.m_mu2=cdip->Lorentz()->Scale(cs.m_z,cs.m_y,cs.m_kt2,Q2);
  cs.m_mu2=Max(cs.m_mu2,cs.m_mode&1?
	       p_shower->GetSudakov()->ISPT2Min():
	       p_shower->GetSudakov()->FSPT2Min());
  if (m_kfmode && !(kmode&16))
  cs.m_mu2*=cdip->Coupling()->CplFac(cs.m_mu2);
  cs.m_cpl=cdip->PureQCD()?0:-1;
  if (!(kmode&2)) return;
  else {
  double eta=1.0;
  if (cs.m_mode==1) eta=GetX(i,cdip)*cs.m_z;
  else if (cs.m_mode==2) eta=GetX(k,cdip)*(1.0-cs.m_y);
  else if (cs.m_mode==3) eta=GetX(i,cdip)*cs.m_z;
  cs.m_wk=(*cdip)(cs.m_z,cs.m_y,eta,-1.0,Q2)*Q2/cs.m_kt2;
  if (IsBad(cs.m_wk)) cs.m_wk=0.0;
  SF_Lorentz *lf = cdip->Lorentz();
  SF_Coupling *cf = cdip->Coupling();
  msg_Debugging()<<"Kernel weight (NLO="<<(kmode&16)
		 <<") [m="<<cs.m_mode<<",c="<<cs.m_col<<"] ( x = "<<eta
		 <<" ) "<<Demangle(typeid(*lf).name()).substr(10)
		 <<"|"<<Demangle(typeid(*cf).name()).substr(10)
		 <<" {\n  "<<*i<<"\n  "<<*j<<"\n  "<<*k
		 <<"\n} -> w = "<<cs.m_wk<<"\n";
  }
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,const ATOOLS::Mass_Selector *const ms,
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
    if (mfb2==0.0 || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
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
    kin=0;
  }
  Kin_Args lt;
  if (i>1) {
    if (k>1) lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,2|(kin?4:0));
    else lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,2|(kin?4:0));
    if (lt.m_pk[3]*(k>1?pk[3]:-pk[3])<0.0 ||
	lt.m_pk[0]<0.0) return Vec4D_Vector();
  }
  else {
    if (k>1) {
      lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-pb,2|(kin?4:0));
      if ((kmode&1) && lt.m_mode) lt.m_stat=-1;
    }
    else lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,2|(kin?4:0));
    if (lt.m_pi[3]*pi[3]>0.0 || lt.m_pi[0]<0.0) return Vec4D_Vector();
  }
  if (lt.m_stat<0) return Vec4D_Vector();
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)j) continue;
    if (m==(size_t)i) after[l]=i>1?lt.m_pi:-lt.m_pi;
    else if (m==(size_t)k && sk) after[l]=k>1?lt.m_pk:-lt.m_pk;
    else after[l]=lt.m_lam*ampl.Leg(m)->Mom();
    ++l;
  }
  if (ampl.Next()) ampl.Next()->SetKin(lt.m_mode&1);
  return after;
}

namespace CSSHOWER {

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs)
  {
    return str<<"CS{kt="<<sqrt(cs.m_kt2)<<",z="<<cs.m_z<<",phi="<<cs.m_phi
	      <<",mode="<<cs.m_mode<<",kin="<<cs.m_kin
	      <<",kmode="<<cs.m_kmode<<"}";
  }

}
