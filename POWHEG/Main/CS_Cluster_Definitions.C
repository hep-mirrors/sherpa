#include "POWHEG/Main/CS_Cluster_Definitions.H"

#include "POWHEG/Showers/Shower.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace POWHEG;
using namespace PDF;
using namespace ATOOLS;

const double s_uxeps=1.0e-3;

CS_Cluster_Definitions::CS_Cluster_Definitions
(Shower *const shower,const int kmode):
  p_shower(shower), m_mode(0), m_kmode(kmode), m_amode(0) {}

CParam CS_Cluster_Definitions::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms)
{
  m_mode=m_kmode;
  CS_Parameters cs(KT2(&ampl,ampl.Leg(i),ampl.Leg(j),ampl.Leg(k),mo,ms));
  m_mode=0;
  return CParam(cs.m_kt2,cs.m_ws,cs.m_x,cs.m_mu2,cs.m_kin);
}

CS_Parameters CS_Cluster_Definitions::KT2
(const ATOOLS::Cluster_Amplitude *ampl,
 const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo,
 ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
  if ((i->Id()&3)<(j->Id()&3)) std::swap<const Cluster_Leg*>(i,j);
  p_b=ampl->Leg(i==ampl->Leg(0)?1:0);
  if ((i->Id()&3)==0) {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) return KT2_FF(i,j,k,mo);
      return KT2_FI(i,j,k,mo);
    }
  }
  else {
    if ((j->Id()&3)==0) {
      if ((k->Id()&3)==0) return KT2_IF(i,j,k,mo);
      return KT2_II(i,j,k,mo);
    }
  }
  THROW(fatal_error,"Unknown CS dipole configuration");  
}

double CS_Cluster_Definitions::Lambda
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

double CS_Cluster_Definitions::Phi
(Vec4D pijt,Vec4D pkt,Vec4D pi,const bool ii) const
{
  Vec4D ktt(0.0,cross(Vec3D(pijt),Vec3D(pkt)));
  Poincare cms(pijt+pkt);
  cms.Boost(pijt);
  cms.Boost(pi);
  Poincare zax(pijt,Vec4D::ZVEC);
  if (!ii && ktt.PSpat2()>1.0e-6) {
    zax.Rotate(ktt);
  }
  else {
    msg_Debugging()<<"Set fixed n_perp\n";
    if (IsZero(pijt[1],1.0e-6) &&
	IsZero(pijt[2],1.0e-6)) pijt[1]=pijt[2]=0.0;
    zax=Poincare(pijt,Vec4D::ZVEC);
    ktt=Vec4D(0.0,1.0,1.0,0.0);
  }
  zax.Rotate(pi);
  Poincare xax(ktt,Vec4D::XVEC);
  xax.Rotate(pi);
  return pi.Phi();
}

double CS_Cluster_Definitions::GetX
(const Cluster_Leg *l,Splitting_Function_Base *const sf) const
{
  const Vec4D &p(l->Mom());
  if (p.PPlus()<p.PMinus()) {
    if (sf) sf->Lorentz()->SetBeam(0);
    return -p.PPlus()/rpa.gen.PBeam(0).PPlus();
  }
  if (sf) sf->Lorentz()->SetBeam(1);
  return -p.PMinus()/rpa.gen.PBeam(1).PMinus();
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
  const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFMap());
  if (cs.m_mode==2) cmap=&p_shower->GetSudakov()->FIMap();
  else if (cs.m_mode==1) cmap=&p_shower->GetSudakov()->IFMap();
  else if (cs.m_mode==3) cmap=&p_shower->GetSudakov()->IIMap();
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
                   <<ProperFlav(k->Flav())<<"("<<k->Flav()<<")].\n";
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
  cs.p_sf=cdip;
  cs.m_mu2=cs.m_kt2*cdip->Coupling()->CplFac(cs.m_kt2);
  cs.m_idi=i->Id();
  cs.m_idj=j->Id();
  cs.m_idk=k->Id();
  if (!(m_mode&1)) return;
  p_shower->SetMS(p_ms);
  cdip->SetFlavourSpec(fls);
  double Q2=dabs((i->Mom()+j->Mom()+k->Mom()).Abs2());
  Parton *pb(NULL);
  double scale=cs.m_kt2, eta=1.0;
  if (cs.m_mode==1) eta=GetX(i,cdip)*cs.m_z;
  else if (cs.m_mode==2) eta=GetX(k,cdip)*(1.0-cs.m_y);
  else if (cs.m_mode==3) {
    eta=GetX(i,cdip)*cs.m_z;
    pb = new Parton(p_b->Flav().Bar(),-cs.m_pbt,pst::IS);
    pb->SetXbj(GetX(p_b,NULL)*cs.m_x/cs.m_z);
    cdip->SetSpec(pb);
  }
  cs.m_wk=(*cdip)(cs.m_z,cs.m_y,eta,scale,Q2);
  if (cs.m_wk<=0.0 || IsBad(cs.m_wk) || 
      (m_amode==1 && !cdip->On()))
    cs.m_wk=sqrt(std::numeric_limits<double>::min());
  cs.m_ws=cs.m_kt2/cs.m_wk;
  msg_Debugging()<<"Kernel weight ["<<cs.m_mode<<"] ( x = "<<eta
		 <<" ) {\n  "<<*i<<"\n  "<<*j<<"\n  "<<*k
		 <<"\n} -> w = "<<cs.m_wk<<" ("<<cs.m_ws<<")\n";
  if (pb) delete pb;
}

CS_Parameters CS_Cluster_Definitions::KT2_FF
(const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo) 
{
  double pipj = i->Mom()*j->Mom();
  double pipk = i->Mom()*k->Mom();
  double pjpk = j->Mom()*k->Mom();
  
  double yijk = pipj/(pipj+pipk+pjpk);
  double zi   = pipk/(pipk+pjpk);

  Vec4D  Q   = i->Mom()+j->Mom()+k->Mom();
  double Q2  = Q*Q;
  
  double mi2 = sqr(p_ms->Mass(i->Flav()));
  double mj2 = sqr(p_ms->Mass(j->Flav()));
  double mk2 = sqr(p_ms->Mass(k->Flav()));
  double mij2 = sqr(p_ms->Mass(mo));
  double pipj2 = (i->Mom()+j->Mom()).Abs2();
  
  double kt2 = (Q2-mi2-mj2-mk2)*yijk*zi*(1.-zi)-(1.-zi)*(1.-zi)*mi2 - zi*zi*mj2;

  double lrat = Lambda(Q2,mij2,mk2)/Lambda(Q2,pipj2,mk2);
  Vec4D pkt(sqrt(lrat)*(k->Mom()-(Q*k->Mom())/Q2*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q);

  if (lrat<0.0 || pkt[0]<0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,0.0,0);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }

  CS_Parameters cs(kt2,zi,yijk,Phi(Q-pkt,pkt,i->Mom()),1.0,Q2,0);
  KernelWeight(i,j,k,mo,cs);
  return cs;
}

CS_Parameters CS_Cluster_Definitions::KT2_FI
(const ATOOLS::Cluster_Leg *i,const ATOOLS::Cluster_Leg *j,
 const ATOOLS::Cluster_Leg *a,const ATOOLS::Flavour &mo) 
{
  // assume pa is outgoing
  double pipj = i->Mom()*j->Mom();
  double pipa = i->Mom()*a->Mom();
  double pjpa = j->Mom()*a->Mom();
  
  double mi2  = sqr(p_ms->Mass(i->Flav()));
  double mj2  = sqr(p_ms->Mass(j->Flav()));
  double mij2 = sqr(p_ms->Mass(mo));
  double ma2  = sqr(p_ms->Mass(a->Flav()));
  
  double xija = (pipa+pjpa+pipj)/(pipa+pjpa);
  double zi   = pipa/(pipa+pjpa);

  double kt2 = -2.*(pipa+pjpa)*(1.-xija)*zi*(1.0-zi)-sqr(1.0-zi)*mi2-zi*zi*mj2;
  if (kt2<0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,0.0,2);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  
  Vec4D Q(a->Mom()+i->Mom()+j->Mom());

  double sij((i->Mom()+j->Mom()).Abs2()), Q2(Q.Abs2());
  double lrat((sqr(Q2-mij2-ma2)-4.0*mij2*ma2)/
	      (sqr(Q2-sij-ma2)-4.0*sij*ma2));
  Vec4D pat(sqrt(lrat)*(a->Mom()-(Q*a->Mom()/Q2)*Q)+(Q2+ma2-mij2)/(2.*Q2)*Q);
  if (lrat<0.0 || pat[0]>0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,0.0,2);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }

  CS_Parameters cs(kt2,zi,1.0-xija,Phi(Q-pat,-pat,i->Mom()),xija,Q2,2);
  KernelWeight(i,j,a,mo,cs);
  return cs;
}

CS_Parameters CS_Cluster_Definitions::KT2_IF
(const ATOOLS::Cluster_Leg *a,const ATOOLS::Cluster_Leg *i,
 const ATOOLS::Cluster_Leg *k,const ATOOLS::Flavour &mo) 
{
  //assume pa is outgoing
  double pipa = i->Mom()*a->Mom();
  double pkpa = k->Mom()*a->Mom();
  double pipk = i->Mom()*k->Mom();
  
  double xika = (pipa+pkpa+pipk)/(pipa+pkpa);
  double ui   = pipa/(pipa+pkpa);
  
  double mi2  = sqr(p_ms->Mass(i->Flav()));
  double mk2  = sqr(p_ms->Mass(k->Flav()));
  double ma2  = sqr(p_ms->Mass(a->Flav()));
  double mai2 = sqr(p_ms->Mass(mo));

  double kt2  = -2.*(pipa+pkpa)*(1.-xika)*ui-mi2-sqr(1.0-xika)*ma2; 
  if (kt2<0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,0.0,1);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  
  Vec4D Q(a->Mom()+i->Mom()+k->Mom()), pi(i->Mom());
  double Q2=Q.Abs2();

  if ((p_shower->KinScheme()==1 ||
       IsZero(xika-ui,s_uxeps)) && ma2==mai2) {
  double sik((i->Mom()+k->Mom()).Abs2()), Q2(Q.Abs2());
  double lrat((sqr(Q2-mk2-ma2)-4.0*mk2*ma2)/
	      (sqr(Q2-sik-ma2)-4.0*sik*ma2));
  Vec4D pat(sqrt(lrat)*(a->Mom()-(Q*a->Mom()/Q2)*Q)+(Q2+ma2-mk2)/(2.*Q2)*Q);
  if (lrat<0.0 || pat[0]>0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,1,1);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  CS_Parameters cs(kt2,xika,ui,Phi(Q-pat,-pat,i->Mom()),xika,Q2,1,1);
  KernelWeight(a,i,k,mo,cs);
  return cs;
  }
  double lrat=Lambda(Q2,mai2,mk2)/Lambda(Q2,(a->Mom()+i->Mom()).Abs2(),mk2);
  Vec4D pkt=sqrt(lrat)*(k->Mom()-(Q*k->Mom()/Q2)*Q)+(Q2+mk2-mai2)/(2.*Q2)*Q;
  if (lrat<0.0 || pkt[0]<0.0) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,0.0,1);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  ZAlign lt(pkt-Q,-p_b->Mom(),mai2,p_ms->Mass2(p_b->Flav()));
  lt.Align(pi);
  lt.Align(pkt);
  CS_Parameters cs(kt2,xika,ui,Phi(lt.PaNew(),pkt,pi),xika,Q2,1);
  if (IsZero(xika-ui,s_uxeps) && ma2==mai2)
    THROW(fatal_error,"Invalid parton configuration");
  KernelWeight(a,i,k,mo,cs);
  return cs;
}

CS_Parameters CS_Cluster_Definitions::KT2_II
(const ATOOLS::Cluster_Leg *a,const ATOOLS::Cluster_Leg *i,
 const ATOOLS::Cluster_Leg *b,const ATOOLS::Flavour &mo) 
{
  //assume pa & pb are outgoing
  double papb = a->Mom()*b->Mom();
  double pipa = i->Mom()*a->Mom();
  double pipb = i->Mom()*b->Mom();
  
  double xiab = (papb+pipa+pipb)/papb;
  double vi   = -pipa/papb;

  double mi2  = sqr(p_ms->Mass(i->Flav()));
  double ma2  = sqr(p_ms->Mass(a->Flav()));
  double mb2  = sqr(p_ms->Mass(b->Flav()));
  double mai2 = sqr(p_ms->Mass(mo));

  double kt2   = 2.*papb*vi*(1.-xiab-vi)-mi2-sqr(1.-xiab-vi)*ma2;

  Vec4D Q(a->Mom()+i->Mom()+b->Mom());
  double Q2=Q.Abs2();

  if (p_shower->KinScheme()==1) {
    double ttau  = Q2-mai2-mb2, tau = Q2-ma2-mi2-mb2;
    double xiiab = xiab*(ttau+sqrt(ttau*ttau-4.*mai2*mb2))/
      (tau+sqrt(tau*tau-4.*ma2*mb2*xiab*xiab));
    double gam   = papb+sqrt(papb*papb-ma2*mb2);
    Vec4D pait(xiiab*(1.0-mai2*mb2/sqr(gam*xiiab))/(1.0-ma2*mb2/sqr(gam))
	       *(-a->Mom()+ma2/gam*b->Mom())-mai2/(xiiab*gam)*b->Mom());
    if (xiab<0.0) {
      CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,3);
      cs.m_wk=cs.m_ws=-1.0;
      return cs;
    }
    CS_Parameters cs(kt2,xiab+vi,vi/xiab,Phi(pait,-b->Mom(),i->Mom(),true),xiab,Q2,3,1);
    KernelWeight(a,i,b,mo,cs);
    return cs;
  }

  double lrat=Lambda(Q2,mai2,mb2)/Lambda(Q2,(a->Mom()+i->Mom()).Abs2(),mb2);
  Vec4D pbt=sqrt(lrat)*(b->Mom()-(Q*b->Mom()/Q2)*Q)+(Q2+mb2-mai2)/(2.*Q2)*Q;
  if (lrat<0.0 || pbt[0]>0.0 || xiab<0.0 ||
      IsZero(sqr(pbt[0])/Q2) || IsZero(sqr(Q[0]-pbt[0])/Q2)) {
    CS_Parameters cs(sqrt(std::numeric_limits<double>::max()),1.0,1.0,0.0,0.0,3);
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  Vec4D pait=-Q+pbt, pb=-b->Mom(), pi=i->Mom();
  ZAlign lt(pait,pb,mai2,mb2);
  pait=lt.PaNew();
  lt.Align(pi);

  CS_Parameters cs(kt2,xiab+vi,vi/xiab,Phi(pait,pb,pi,true),xiab,Q2,3);
  if (pait[3]*pb[3]>0.0) {
    cs.m_wk=cs.m_ws=-1.0;
    return cs;
  }
  cs.m_pbt=pbt;
  lt.Align(cs.m_pbt);
  KernelWeight(a,i,b,mo,cs);
  return cs;
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo,ATOOLS::Mass_Selector *const ms,
 const int kin)
{
  p_ms=ms;
  if (i>j) std::swap<int>(i,j);
  if (i>1 && j>1 && k>1) return Combine_FF(ampl,i,j,k,mo);
  if (i>1 && j>1 && k<2) return Combine_FI(ampl,i,j,k,mo);
  if (i<2 && j>1 && k>1) return Combine_IF(ampl,i,j,k,mo,kin);
  if (i<2 && j>1 && k<2) return Combine_II(ampl,i,j,k,mo,kin);
  
  std::cout<<" asked for combine of : "<<i<<" "<<j<<" "<<k<<std::endl;  
  THROW(fatal_error,"Unknown CS dipole configuration"); 
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine_FF
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const ATOOLS::Flavour &mo)
{
  if (i>j) std::swap<int>(i,j);
  Vec4D_Vector after(ampl.Legs().size()-1);
 
  //old momenta
  Vec4D pi = ampl.Leg(i)->Mom();
  Vec4D pj = ampl.Leg(j)->Mom();
  Vec4D pk = ampl.Leg(k)->Mom();
  Vec4D Q  = pi+pj+pk;
  
  //masses etc.
  double mij2 = sqr(p_ms->Mass(mo));
  double mk2  = sqr(p_ms->Mass(ampl.Leg(k)->Flav()));
  double Q2   = Q*Q;
  
  //new momenta
  double lrat = Lambda(Q2,mij2,mk2)/Lambda(Q2,(pi+pj)*(pi+pj),mk2);
  Vec4D pkt  = sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+ (Q2+mk2-mij2)/(2.*Q2)*Q;
  if (lrat<0.0 || pkt[0]<0.0) return Vec4D_Vector();
  Vec4D pijt = Q-pkt; 

  //setting the new momenta
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)j) continue;
    after[l]=ampl.Leg(m)->Mom();
    if (m==(size_t)i) after[l]=pijt;
    if (m==(size_t)k) after[l]=pkt;
    ++l;
  }
  return after;
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine_FI
(const Cluster_Amplitude &ampl,int i,int j,int a,
 const ATOOLS::Flavour &mo)
{
  if (i>j) std::swap<int>(i,j);
  Vec4D_Vector after(ampl.Legs().size()-1);
  
  //old momenta (assume pa is outgoing)
  Vec4D pi = ampl.Leg(i)->Mom();
  Vec4D pj = ampl.Leg(j)->Mom();
  Vec4D pa = ampl.Leg(a)->Mom();
  
  //masses etc
  double mij2 = p_ms->Mass2(mo);
  double ma2  = p_ms->Mass2(ampl.Leg(a)->Flav());
  double mb2  = p_ms->Mass2(ampl.Leg(1-a)->Flav());

  Vec4D Q(pa+pi+pj);

  double sij((pi+pj).Abs2()), Q2(Q.Abs2());
  double lrat((sqr(Q2-mij2-ma2)-4.0*mij2*ma2)/
	      (sqr(Q2-sij-ma2)-4.0*sij*ma2));
  Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+ma2-mij2)/(2.*Q2)*Q);
  if (lrat<0.0 || pat[0]>0.0) return Vec4D_Vector();
  Vec4D pijt(Q-pat);

  Vec4D pb=ampl.Leg(1-a)->Mom();
  if (pat[3]*pb[3]>0.0) return Vec4D_Vector();
  ZAlign lt(-pat,-pb,ma2,mb2);
  Vec4D pan(-lt.PaNew());
  if (pan[0]>0.0 || IsZero(pan[0],1.0e-6) ||
      lt.Status()<0) return Vec4D_Vector();
  if (-pan[0]>rpa.gen.PBeam(a)[0]) return Vec4D_Vector();
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)j) continue;
    if (m==(size_t)a) after[l]=pan;
    else if (m==(size_t)1-a) after[l]=pb;
    else {
      after[l]=ampl.Leg(m)->Mom();
      if (m==(size_t)i) after[l]=pijt;
      lt.Align(after[l]);
    }
    ++l;
  }
  return after;
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine_IF
(const Cluster_Amplitude &ampl,int a,int i,int k,
 const ATOOLS::Flavour &mo,const int kin)
{
  Vec4D_Vector after(ampl.Legs().size()-1);
  
  //old momenta (assume pa is outgoing)
  Vec4D pa = ampl.Leg(a)->Mom();
  Vec4D pi = ampl.Leg(i)->Mom();
  Vec4D pk = ampl.Leg(k)->Mom();
  
  double pipa = pi*pa;
  double pkpa = pk*pa;
  double pipk = pi*pk;
  
  double xika = (pipa+pkpa+pipk)/(pipa+pkpa);
  double ui   = pipa/(pipa+pkpa);
  
  double mi2  = p_ms->Mass2(ampl.Leg(i)->Flav());
  double mk2  = p_ms->Mass2(ampl.Leg(k)->Flav());
  double ma2  = p_ms->Mass2(ampl.Leg(a)->Flav());
  double mai2 = p_ms->Mass2(mo);
  double mb2  = p_ms->Mass2(ampl.Leg(1-a)->Flav());

  Vec4D Q(pa+pi+pk);
  double Q2=Q.Abs2();

  if (kin==1) {
  double kt2  = -2.*(pipa+pkpa)*(1.-xika)*ui-mi2-sqr(1.0-xika)*ma2; 
  if (kt2<0.0) return Vec4D_Vector();

  double sik((pi+pk).Abs2()), Q2(Q.Abs2());
  double lrat((sqr(Q2-mk2-ma2)-4.0*mk2*ma2)/
	      (sqr(Q2-sik-ma2)-4.0*sik*ma2));
  Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+ma2-mk2)/(2.*Q2)*Q);
  if (lrat<0.0 || pat[0]>0.0) return Vec4D_Vector();
  Vec4D pikt(Q-pat);

  Vec4D pb=ampl.Leg(1-a)->Mom();
  if (pat[3]*pb[3]>0.0) return Vec4D_Vector();
  ZAlign lt(-pat,-pb,ma2,mb2);
  Vec4D pan(-lt.PaNew());
  if (pan[0]>0.0 || IsZero(pan[0],1.0e-6) ||
      lt.Status()<0) return Vec4D_Vector();
  if (-pan[0]>rpa.gen.PBeam(a)[0]) return Vec4D_Vector();
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)i) continue;
    if (m==(size_t)a) after[l]=pan;
    else if (m==(size_t)1-a) after[l]=pb;
    else {
      after[l]=ampl.Leg(m)->Mom();
      if (m==(size_t)k) after[l]=pikt;
      lt.Align(after[l]);
    }
    ++l;
  }
  return after;
  }
  double lrat=Lambda(Q2,mai2,mk2)/Lambda(Q2,(pa+pi)*(pa+pi),mk2);
  Vec4D pikt=sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mai2)/(2.*Q2)*Q;
  if (lrat<0.0 || pikt[0]<0.0) return Vec4D_Vector();
  Vec4D pb(ampl.Leg(1-a)->Mom());
  ZAlign lt(pikt-Q,-pb,mai2,mb2);
  Vec4D pan(-lt.PaNew());
  if (pan[0]>0.0 || IsZero(pan[0],1.0e-6) ||
      lt.Status()<0) return Vec4D_Vector();
  if (-pan[0]>rpa.gen.PBeam(a)[0]) return Vec4D_Vector();
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)i) continue;
    if (m==(size_t)a) after[l]=pan;
    else if (m==(size_t)1-a) after[l]=pb;
    else {
      after[l]=ampl.Leg(m)->Mom();
      if (m==(size_t)k) after[l]=pikt;
      lt.Align(after[l]);
    }
    ++l;
  }
  return after;
}

ATOOLS::Vec4D_Vector  CS_Cluster_Definitions::Combine_II
(const Cluster_Amplitude &ampl,int a,int i,int b,
 const ATOOLS::Flavour &mo,const int kin)
{
  Vec4D_Vector after(ampl.Legs().size()-1);
  
  //old momenta (assume pa & pb are outgoing)
  Vec4D pa = ampl.Leg(a)->Mom();
  Vec4D pi = ampl.Leg(i)->Mom();
  Vec4D pb = ampl.Leg(b)->Mom();
  
  Vec4D Q(pa+pi+pb), pai(pa+pi);

  double mb2  = sqr(p_ms->Mass(ampl.Leg(b)->Flav()));
  double mai2 = sqr(p_ms->Mass(mo));
  double sai  = pai.Abs2(), Q2 = Q.Abs2();

  if (kin==1) {
  double papb = pa*pb, pipa = pi*pa, pipb = pi*pb;
  double xiab = (papb+pipa+pipb)/papb;
  if (xiab<0.0) return Vec4D_Vector();

  double mi2  = sqr(p_ms->Mass(ampl.Leg(i)->Flav()));
  double ma2  = sqr(p_ms->Mass(ampl.Leg(a)->Flav()));
  double Q2   = (pa+pi+pb).Abs2();

  double ttau  = Q2-mai2-mb2, tau = Q2-ma2-mi2-mb2;
  double xiiab = xiab*(ttau+sqrt(ttau*ttau-4.0*mai2*mb2))
    /(tau+sqrt(tau*tau-4.0*ma2*mb2*xiab*xiab));
  if (ttau*ttau<4.0*mai2*mb2 || ttau<=0.0 ||
      tau*tau<4.0*ma2*mb2*xiab*xiab || tau<=0.0) return Vec4D_Vector();
  double gam   = papb+sqrt(papb*papb-ma2*mb2);
  if (IsZero(xiiab,1.0e-6) || IsZero(gam,1.0e-6)) return Vec4D_Vector();

  Vec4D pait = xiiab
    *(1.0-mai2*mb2/sqr(gam*xiiab))/(1.0-ma2*mb2/sqr(gam))
    *(pa-ma2/gam*pb)+mai2/(xiiab*gam)*pb;
  if (pait[3]*pb[3]>0.0) return Vec4D_Vector();
  Vec4D pan = pait;

  Poincare cmso(-pa-pb-pi), cmsn(-pait-pb);
  cmso.Boost(pait);
  Poincare zrot(pait,-Sign(pb[3])*Vec4D::ZVEC);
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)i) continue;
    if (m==(size_t)a) after[l]=pan;
    else if (m==(size_t)1-a) after[l]=pb;
    else {
      after[l]=ampl.Leg(m)->Mom();
      cmso.Boost(after[l]);
      zrot.Rotate(after[l]);
      cmsn.BoostBack(after[l]);
    }
    ++l;
  }
  return after;
  }

  ZAlign lt(-pai,-pb,sai,mb2);
  pai=-lt.PaNew();
  Q=pai+pb;
  if (pai[3]*pb[3]>0.0) return Vec4D_Vector();

  double lrat=Lambda(Q2,mai2,mb2)/Lambda(Q2,sai,mb2);
  Vec4D pbt=sqrt(lrat)*(pb-(Q*pb/Q2)*Q)+(Q2+mb2-mai2)/(2.*Q2)*Q, pait=Q-pbt;
  double xiab = (pa*pb+pi*pa+pi*pb)/(pa*pb);
  if (lrat<0.0 || pbt[0]>0.0 || 
      IsZero(sqr(pbt[0])/Q2) || xiab<0.0) return Vec4D_Vector();

  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)i) continue;
    if (m==(size_t)a) after[l]=pait;
    else if (m==(size_t)1-a) after[l]=pbt;
    else {
      after[l]=ampl.Leg(m)->Mom();
      lt.Align(after[l]);
    }
    ++l;
  }
  return after;
}

namespace POWHEG {

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs)
  {
    return str<<"CS{kt="<<sqrt(cs.m_kt2)<<",z="<<cs.m_z<<",phi="<<cs.m_phi
	      <<",mode="<<cs.m_mode<<",kin="<<cs.m_kin<<"}";
  }

}
