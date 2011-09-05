#include "PHASIC++/Channels/VHAAG_res.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;

const double a1_s  =1.5;
const double s_s   =.3;
const double a1_sp =1.;
const double s_sp1 =1.;
const double s_sp2 =.5;
const double a1_i  =1.;
const double a2_i  =1.;
const double s_i   =1.;
const double a1_iF =1.;
const double a2_iF =.0;
const double s_xx  =.5;

VHAAG_res::VHAAG_res(int _nin,int _nout,int pn,VHAAG_res* ovl)
{
  nin=_nin; nout=_nout;n_ap=nin+nout-1;
  if (n_ap<5) {
    msg_Error()<<"Minimum number of final state particles for VHAAG_res integrator is 4!"<<std::endl;
    abort();
  }
  if (!ovl) {
    Data_Reader dr(" ",";","!","=");
    dr.AddWordSeparator("\t");
    dr.SetInputPath(rpa->GetPath());
    dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
    m_rkf = dr.GetValue<int>("VHAAG_RES_KF",24);
    n_d1  = dr.GetValue<int>("VHAAG_RES_D1",2);
    n_d2  = dr.GetValue<int>("VHAAG_RES_D2",3);
    msg_Out()<<" Initialized VHAAG with "<<Flavour((kf_code)(m_rkf))<<"-resonance."<<std::endl;
  }
  else {
    m_rkf = ovl->RKF();
    n_d1  = ovl->RD1();
    n_d2  = ovl->RD2();
  }

  Permutation pp(n_ap-2);
  int bp=pn%2;
  pn/=2;
  int* tp=pp.Get(pn);
  int i;
  for (i=0;i<n_ap-2;i++) tp[i]++;
  for (i=0;i<n_ap-2;i++) if (tp[i]>=n_d1) tp[i]+=2;
  std::vector<int> perm(n_ap);
  perm[0] = 0;
  if (bp==0) {
    for (i=1;tp[i-1]!=1;i++) perm[i] = tp[i-1];
    n_b=i;perm[i]=n_d1;
    for (;i<n_ap-1;i++) perm[i+1] = tp[i-1];
  }
  else {
    for (i=1;i<n_ap-1;i++) perm[i] = tp[i-1];
    n_b=n_ap-1;perm[n_b]=n_d1;
  }
  Initialize(perm,ovl);
}

void VHAAG_res::Initialize(std::vector<int> perm, VHAAG_res* ovl) 
{
  p_perm = new int[n_ap];
  p_perm[0] = 0;
  p_s = new double[n_ap];
  for (int i=0;i<n_ap;i++) p_s[i]=0.;
  msg_Tracking()<<"Init VHAAG_res: 0";
  name = "VHAAG_res";
  char hlp[4];
  for (int i=1;i<n_ap;i++) {
    p_perm[i] = perm[i];
    if (perm[i]==1) n_p1=i;
    name+= "_";
    sprintf(hlp,"%i",p_perm[i]);
    name+= std::string(hlp);
    msg_Tracking()<<" "<<p_perm[i];
  }

  rannum = 3*nout-4;;
  rans  = new double[rannum];
  m_q = new Vec4D[n_ap];
  m_ownvegas = false;
  if (ovl) p_sharedvegaslist = ovl->GetSharedVegasList();
  else p_sharedvegaslist = NULL;
  if (p_sharedvegaslist==NULL) {
    p_sharedvegaslist = new std::map<int,Vegas*>;
  }

  m_type=Min(n_p1-1,n_ap-(n_p1+1))+1;
  if (n_p1-1>=n_ap-(n_p1+1) && n_b>n_p1) m_type=-m_type;
  if (n_p1-1<n_ap-(n_p1+1) && n_b<n_p1) m_type=-m_type;
  msg_Tracking()<<" n_p1="<<n_p1<<" type="<<m_type<<std::endl;

  if (1) {
    if (p_sharedvegaslist->find(m_type)==p_sharedvegaslist->end()) {
      (*p_sharedvegaslist)[m_type] = new Vegas(rannum,100,Name());
//       (*p_sharedvegaslist)[m_type]->SetAutoOptimize(n_ap*10);//Min(nout,5)*15);
      if (0) {
	if (abs(m_type)<3) {
	  for (int i=0;i<=n_ap-5;i++) (*p_sharedvegaslist)[m_type]->ConstChannel(2+3*i);
	} 
	else {
	  for (int i=0;i<=abs(m_type)-3;i++) (*p_sharedvegaslist)[m_type]->ConstChannel(3+3*i);
	  (*p_sharedvegaslist)[m_type]->ConstChannel(3*(abs(m_type)-2)+2);
	  for (int i=0;i<n_ap-abs(m_type)-3;i++) 
	    (*p_sharedvegaslist)[m_type]->ConstChannel(3*abs(m_type)+3*i-1);
	} 
	(*p_sharedvegaslist)[m_type]->ConstChannel(rannum-4);
      } 
      m_ownvegas = true;
    } 
    p_vegas = (*p_sharedvegaslist)[m_type];
  } 
  else  {
    m_ownvegas = true;
    p_vegas = new Vegas(rannum,100,Name());
    p_vegas->SetAutoOptimize(50);
    delete p_sharedvegaslist;
    p_sharedvegaslist=0;
    if (0) {
      if (abs(m_type)<3) {
	for (int i=0;i<=n_ap-5;i++) p_vegas->ConstChannel(2+3*i);
      } 
      else {
	for (int i=0;i<=abs(m_type)-3;i++) p_vegas->ConstChannel(3+3*i);
	p_vegas->ConstChannel(3*(abs(m_type)-2)+2);
	for (int i=0;i<n_ap-abs(m_type)-3;i++) p_vegas->ConstChannel(3*abs(m_type)+3*i-1);
      } 
      p_vegas->ConstChannel(rannum-4);
    } 
  } 

  m_s0=-1.;
}

VHAAG_res::~VHAAG_res()
{
  delete[] p_perm;
  delete[] p_s;
  delete[] m_q;
  if (m_ownvegas) {
    delete p_vegas;
    if (p_sharedvegaslist) p_sharedvegaslist->erase(m_type);
  }
  if (p_sharedvegaslist) if (p_sharedvegaslist->empty()) {
    delete p_sharedvegaslist;
  }
}

void VHAAG_res::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}

double VHAAG_res::PiFunc(double a1,double a2,
		     double s1b,double s2b,double c)
{
  return 4.*(1.-c*c)*((1.-a2+s2b-s1b)*a2-s2b)
    -sqr((1.-2.*a1+s1b-s2b)+(1.-2.*a2-s1b+s2b)*c);
}

void VHAAG_res::Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n1,int n2,double *ran)
{
  double s = Q.Abs2();
  double s1min = 0.;
  double s2min = 0.;
  if (n1<n2) {
    for (int i=n1+1;i<n2;i++) s1min+=p_s[i];
    for (int i=n2+1;i<n_ap;i++) s2min+=p_s[i];
  }
  else {
    for (int i=n1+1;i<n_ap;i++) s1min+=p_s[i];
    for (int i=n2+1;i<n1;i++) s2min+=p_s[i];
  }
  double s1max = sqr(sqrt(s)-sqrt(s2min));
  double exp = s_sp1; if (s1min==0.) exp = s_xx;
  double s1 = CE.MasslessPropMomenta(exp,s1min,s1max,ran[0]);

  exp = s_sp2; if (s2min==0.) exp = s_xx;
  double s2max = sqr(sqrt(s)-sqrt(s1));
  double s2 = CE.MasslessPropMomenta(exp,s2min,s2max,ran[1]);

  double pb0 =0.5*(s+s1-s2)/s;
  double pb  =sqrt(pb0*pb0-s1/s);

  double a1min = pb0-pb;
  double a1max = pb0+pb;
  exp = a1_sp; if (a1min==0.) exp = s_xx;
  double a1 = CE.MasslessPropMomenta(exp,a1min,a1max,ran[2]);

  double phi=2.*M_PI*ran[3];

  ConstructMomenta(a1,phi,s1,s2,s,q1,p1,p2);
//     std::cout<<"Split generated: "<<p1<<p2<<std::endl;
//     std::cout<<" s1: "<<s1min<<" < "<<s1<<" < "<<s1max<<std::endl;
//     std::cout<<" s2: "<<s2min<<" < "<<s2<<" < "<<s2max<<std::endl;
//     std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<" "<<a1*(q1*q2)<<std::endl;
}

void VHAAG_res::Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n1,int n2,double *ran)
{
  double s = Q.Abs2();
  double smin = 0.;
  double s1 = p_s[n1];
  for (int i=n2;i<n2+n_ap-3;i++) smin+=p_s[i];
  double smax = sqr(sqrt(s)-sqrt(s1));
  double exp = s_s; if (smin==0.) exp = s_xx;
  double s2 = CE.MasslessPropMomenta(exp,smin,smax,ran[0]);

  double pb0 = 0.5*(s+s1-s2)/s;
  double pb  =sqrt(pb0*pb0-s1/s);
  double a1min = pb0-pb;
  double a1max = pb0+pb;
  double a1 = CE.MasslessPropMomenta(s_xx,a1min,a1max,ran[1]);

  double phi=2.*M_PI*ran[2];

  ConstructMomenta(a1,phi,s1,s2,s,q1,p1,p2);
}

void VHAAG_res::SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double smin,double *ran)
{
  Poincare qb(Q);
  qb.Boost(q1);
  double s = Q.Abs2();
  double smax = s;
  double exp = s_i; if (smin==0.) exp = s_xx;
  double s2 = CE.MasslessPropMomenta(exp,smin,smax,ran[0]);

  double a1min = 0.;
  double a1max = 1.-s2/s;
  exp = a1_i; if (a1min==0.) exp = s_xx;
  double a1 = CE.MasslessPropMomenta(exp,a1min,a1max,ran[1]);

    double phi=2.*M_PI*ran[2];
    Vec4D qq(1.,0.,0.,1.);
    ConstructMomenta(a1,phi,0.,s2,s,qq,p1,p2);
    Poincare rot(qq,q1);
    rot.Rotate(p1);
    rot.Rotate(p2);

  qb.BoostBack(p1);
  qb.BoostBack(p2);
}

void VHAAG_res::SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double s2,double *ran)
{
  Poincare qb(Q);
  qb.Boost(q1);
  double s = Q.Abs2();

  double a1min = 0.;
  double a1max = 1.-s2/s;
  double a1 = CE.MasslessPropMomenta(s_xx,a1min,a1max,ran[0]);

    double phi=2.*M_PI*ran[1];
    Vec4D qq(1.,0.,0.,1.);
    ConstructMomenta(a1,phi,0.,s2,s,qq,p1,p2);
    Poincare rot(qq,q1);
    rot.Rotate(p1);
    rot.Rotate(p2);

  qb.BoostBack(p1);
  qb.BoostBack(p2);
}

void VHAAG_res::SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double s2,double *ran)
{
  double s = Q.Abs2();

  double a1min = 0.; 
  double a1max = 1.-s2/s; 
  double a1 = CE.AntennaMomenta(a1min,a1max,ran[0]);

  double phi=2.*M_PI*ran[1];

  ConstructMomenta(a1,phi,0.,s2,s,q1,p1,p2);
}

void VHAAG_res::ConstructMomenta(double a1,double a2,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)
{
  double ps = 0.25*(sqr(s-s1-s2)-4.*s1*s2)/s;
  Vec3D e1  = Vec3D(q1)/q1[0];
  Vec3D e2  = Vec3D(q2)/q2[0];
  Vec3D ee  = cross(e1,e2);
  ee = (1./ee.Abs())*ee;
  double v  = e1*e2;
  double v1 = sqrt(ps+s1)-a1*sqrt(s);
  double v2 = sqrt(ps+s2)-a2*sqrt(s);
  double a  = (v1+v2*v)/(1-v*v);
  double b  = -(v2+v1*v)/(1-v*v);
  double eps= sqrt(ps-a*a-b*b-2.*a*b*v);
  if (ran->Get()<0.5) eps=-eps;
  Vec3D pv = a*e1+b*e2+eps*ee;

  p1 = Vec4D(sqrt(ps+s1),pv);
  p2 = Vec4D(sqrt(ps+s2),-1.*pv);
}

void VHAAG_res::ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)
{
  double ps = 0.25*(sqr(s-s1-s2)-4.*s1*s2)/s;
  if (q1.PPerp()!=0.||!IsZero(q1.Abs2()/Max(1.,sqr(q1[0])),1.0e-6)) {
    msg_Error()<<" Error in"<<std::endl
	       <<"ConstructMomenta(double a1,double phi,double s1,double s2,double s,"<<std::endl
	       <<"                 ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)!"<<std::endl
	       <<" q1 must be in beam direction and massless!   q1="<<q1<<" ("<<q1.Abs2()<<")"<<std::endl;
    abort();
  }
  Vec3D e1  = Vec3D(q1)/q1[0];
  double v1 = sqrt(ps+s1)-a1*sqrt(s);

  double cc = sqrt(ps-v1*v1); 
  Vec3D pv(cc*cos(phi),cc*sin(phi),v1*e1[3]);

  p1 = Vec4D(sqrt(ps+s1),pv);
  p2 = Vec4D(sqrt(ps+s2),-1.*pv);
}
    
double VHAAG_res::SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			  ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n1, int n2,double *ran)
{
  double wt=1.;
  double s = Q.Abs2();
  double s1min = 0.;
  double s2min = 0.;
  if (n1<n2) {
    for (int i=n1+1;i<n2;i++) s1min+=p_s[i];
    for (int i=n2+1;i<n_ap;i++) s2min+=p_s[i];
  }
  else {
    for (int i=n1+1;i<n_ap;i++) s1min+=p_s[i];
    for (int i=n2+1;i<n1;i++) s2min+=p_s[i];
  }
  double s1max = sqr(sqrt(s)-sqrt(s2min));
  double exp = s_sp1; if (s1min==0.) exp = s_xx;
  double s1 = p1.Abs2(); 
  wt*= CE.MasslessPropWeight(exp,s1min,s1max,s1,ran[0]);

  exp = s_sp2; if (s2min==0.) exp = s_xx;
  double s2max = sqr(sqrt(s)-sqrt(s1));
  double s2 = p2.Abs2(); 
  wt*= CE.MasslessPropWeight(exp,s2min,s2max,s2,ran[1]);

  double pb0 =0.5*(s+s1-s2)/s;
  double pb  =sqrt(pb0*pb0-s1/s);

  double a1min = pb0-pb;
  double a1max = pb0+pb;
  exp = a1_sp; if (a1min==0.) exp = s_xx;
  double a1 = (q1*p1)/(q1*Q);
  wt*= CE.MasslessPropWeight(exp,a1min,a1max,a1,ran[2]);

  wt*= 2./M_PI;
  ran[3]=p1.Phi()/(2.*M_PI);
  if(ran[3]<0.) ran[3]+=1.;
  return wt;
}

double VHAAG_res::Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			       ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n1,int n2,double *ran)
{
  double wt=1.;
  double s = Q.Abs2();
  double smin = 0.;
  double s1 = p_s[n1];
  for (int i=n2;i<n2+n_ap-3;i++) smin+=p_s[i];
  double smax = sqr(sqrt(s)-sqrt(s1));
  double exp = s_s; if (smin==0.) exp = s_xx;
  double s2 = p2.Abs2(); 
  wt*= CE.MasslessPropWeight(exp,smin,smax,s2,ran[0]);

  double pb0 = 0.5*(s+s1-s2)/s;
  double pb  =sqrt(pb0*pb0-s1/s);
  double a1min = pb0-pb;
  double a1max = pb0+pb;
  double a1 = q1*p1/(q1*Q);
  wt*= CE.MasslessPropWeight(s_xx,a1min,a1max,a1,ran[1]);

  wt*= 2./M_PI;
  ran[2]=p1.Phi()/(2.*M_PI);
  if(ran[2]<0.) ran[2]+=1.;
  return wt;
}

double VHAAG_res::SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
				ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double smin,double *ran)
{
  double wt=1.;
  Q = p1+p2;
  double s = Q.Abs2();
  double s2 = p2.Abs2();
  double smax = s;
  double exp = s_i; if (smin==0.) exp = s_xx;
  wt*= CE.MasslessPropWeight(exp,smin,smax,s2,ran[0]);

  double a1min = 0.;
  double a1max = 1.-s2/s;
  exp = a1_i; if (a1min==0.) exp = s_xx;
  double a1 = q1*p1/(q1*Q);
  wt*= CE.MasslessPropWeight(exp,a1min,a1max,a1,ran[1]);

  Poincare qb(Q);
  qb.Boost(q1);
  wt*=2./M_PI;

    Vec4D qq(1.,0.,0.,1.);
    qb.Boost(p1);
    Poincare rot(qq,q1);
    rot.RotateBack(p1);
    ran[2]=p1.Phi()/(2.*M_PI);
    if(ran[2]<0.) ran[2]+=1.;

  return wt;
}

double VHAAG_res::SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
				 ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double s2,double *ran)
{
  double wt=1.;
  Q = p1+p2;
  double s=Q.Abs2();
  double a1min = 0.;
  double a1max = 1.-s2/s;
  double a1 = q1*p1/(q1*Q);
  wt*= CE.MasslessPropWeight(s_xx,a1min,a1max,a1,ran[0]);

  Poincare qb(Q);
  qb.Boost(q1);

    wt*=2./M_PI;
    Vec4D qq(1.,0.,0.,1.);
    qb.Boost(p1);
    Poincare rot(qq,q1);
    rot.RotateBack(p1);
    ran[1]=p1.Phi()/(2.*M_PI);
    if(ran[1]<0.) ran[1]+=1.;

  return wt;
}

double VHAAG_res::SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
				      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double s2,double *ran)
{
  double wt=1.;
  double a1min = 0.; 
  double a1max = 1.-s2/Q.Abs2(); 
  double a1 = q1*p1/(q1*Q);
  wt*= CE.AntennaWeight(a1min,a1max,a1,ran[0]);

  wt*=2./M_PI;
  ran[1]=p1.Phi()/(2.*M_PI);
  if(ran[1]<0.) ran[1]+=1.;

  return wt;
}

void VHAAG_res::GenerateBosonMass(ATOOLS::Vec4D *p,double *ran)
{
  double smax=(p[0]+p[1]).Abs2();
  Flavour fl((kf_code)(m_rkf));
  p_s[n_b]=CE.MassivePropMomenta(fl.Mass(),fl.Width(),1,0.,smax,ran[rannum-3]);
}

double VHAAG_res::BosonWeight(ATOOLS::Vec4D *p,double *ran)
{
  double smax=(p[0]+p[1]).Abs2();
  Flavour fl((kf_code)(m_rkf));
  double w=CE.MassivePropWeight(fl.Mass(),fl.Width(),1,0.,smax,m_q[n_b].Abs2(),ran[rannum-3]);
  w*=CE.Isotropic2Weight(p[n_d1],p[n_d2],ran[rannum-2],ran[rannum-1]);
  return w;
}


void VHAAG_res::GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			       ATOOLS::Vec4D* q,double* s,int n,double *ran)
{
  Vec4D r1 = q1;
  Vec4D rQ = Q;
  double ms(0.);
  for (int i=1;i<n;i++) ms+=s[i]; 
  for (int i=n;i>2;i--) {
    SingleSplit(r1,rQ,q[n-i],rQ,ms,ran);
    ran+=3;
    r1 = q[n-i];
    ms-= s[n-i+1];
  }
  SingleSplitF(r1,rQ,q[n-2],q[n-1],ms,ran);
}

double VHAAG_res::BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
			 ATOOLS::Vec4D* q,double* s,int n,double *ran)
{
  double wt=1.;
  ran+= 3*(n-2);
  double ms(s[n-1]);
  wt*=SingleSplitFWeight(q[n-3],Q,q[n-2],q[n-1],ms,ran);

  for (int i=3;i<=n;i++) {
    ran-= 3;
    ms+= s[n-i+1];
    wt*=SingleSplitWeight(q[n-i-1],Q,q[n-i],Q,ms,ran);
  }
  return wt;
}

void VHAAG_res::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  CalculateS0(cuts);
  double wt=1.;

  if (n_ap==4) {
    Vec4D Q(p[0]+p[1]);
    wt=SingleSplitF0Weight(p[0],Q,p[2],p[3],p_s[3],rans);  
    weight = p_vegas->GenerateWeight(rans)/wt/pow(2.*M_PI,2);
    return;
  }

  for (int i=0;i<n_ap;i++) m_q[i]=p[p_perm[i]];
  m_q[n_b]=p[n_d1]+p[n_d2];
  p_s[n_b]=m_q[n_b].Abs2();
  Vec4D Q(m_q[0]+m_q[n_p1]);
  
  if (n_p1==1){
    Vec4D P;
    wt*=BranchWeight(m_q[2],P,&(m_q[3]),&(p_s[3]),n_ap-3,rans+3);
    wt*=Split0Weight(m_q[1],Q,m_q[2],P,2,3,rans);    
  }
  else if (n_p1==n_ap-1){
    Vec4D P;
    wt*=BranchWeight(m_q[1],P,&(m_q[2]),&(p_s[2]),n_ap-3,rans+3);
    wt*=Split0Weight(m_q[0],Q,m_q[1],P,1,2,rans);    
  }
  else if (n_p1==2){
    Vec4D P;
    wt*=BranchWeight(m_q[2],P,&(m_q[3]),&(p_s[3]),n_ap-3,rans+3);
    wt*=Split0Weight(m_q[0],Q,m_q[1],P,1,3,rans);    
  }
  else if (n_p1==n_ap-2){
    Vec4D P;
    wt*=BranchWeight(m_q[0],P,&(m_q[1]),&(p_s[1]),n_ap-3,rans+3);
    wt*=Split0Weight(m_q[n_p1],Q,m_q[n_ap-1],P,n_ap-1,1,rans);    
  }
  else if (n_p1<=(n_ap-1)/2) {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[0],Q1,&(m_q[1]),&(p_s[1]),n_p1-1,rans+4);
    wt*=BranchWeight(m_q[n_p1],Q2,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,rans+3*(n_p1-1));
    wt*=SplitWeight(m_q[0],Q,Q1,Q2,0,n_p1,rans);
  }
  else {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[n_p1],Q1,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,rans+4);
    wt*=BranchWeight(m_q[0],Q2,&(m_q[1]),&(p_s[1]),n_p1-1,rans+3*(n_ap-n_p1-1));
    wt*=SplitWeight(m_q[n_p1],Q,Q1,Q2,n_p1,0,rans);
  }
  wt*=BosonWeight(p,rans);
  double vw = p_vegas->GenerateWeight(rans);
  weight = vw/wt/pow(2.*M_PI,nout*3.-4.);
}

void VHAAG_res::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *ran)
{
  CalculateS0(cuts);
  double *vran = p_vegas->GeneratePoint(ran);
  for(int i=0;i<rannum;i++) rans[i]=vran[i];

  GenerateBosonMass(p,vran);
  if (n_ap==4) {
    Vec4D Q(p[0]+p[1]);
    SingleSplitF0(p[0],Q,p[2],p[3],p_s[3],vran);  
    return;
  }

  m_q[0] = p[0];
  m_q[n_p1] = p[1];

  Vec4D Q(m_q[0]+m_q[n_p1]);
  if (n_p1==1){
    Split0(m_q[1],Q,m_q[2],Q,2,3,vran);    
    GenerateBranch(m_q[2],Q,&(m_q[3]),&(p_s[3]),n_ap-3,vran+3);
  }
  else if (n_p1==n_ap-1){
    Split0(m_q[0],Q,m_q[1],Q,1,2,vran);    
    GenerateBranch(m_q[1],Q,&(m_q[2]),&(p_s[2]),n_ap-3,vran+3);
  }
  else if (n_p1==2){
    Split0(m_q[0],Q,m_q[1],Q,1,3,vran);    
    GenerateBranch(m_q[2],Q,&(m_q[3]),&(p_s[3]),n_ap-3,vran+3);
  }
  else if (n_p1==n_ap-2){
    Split0(m_q[n_p1],Q,m_q[n_ap-1],Q,n_ap-1,1,vran);    
    GenerateBranch(m_q[0],Q,&(m_q[1]),&(p_s[1]),n_ap-3,vran+3);
  }
  else if (n_p1<=(n_ap-1)/2) {
    Vec4D Q1,Q2;
    Split(m_q[0],Q,Q1,Q2,0,n_p1,vran);
    GenerateBranch(m_q[0],Q1,&(m_q[1]),&(p_s[1]),n_p1-1,vran+4);
    GenerateBranch(m_q[n_p1],Q2,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,vran+3*(n_p1-1));
  }
  else {
    Vec4D Q1,Q2;
    Split(m_q[n_p1],Q,Q1,Q2,n_p1,0,vran);
    GenerateBranch(m_q[n_p1],Q1,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,vran+4);
    GenerateBranch(m_q[0],Q2,&(m_q[1]),&(p_s[1]),n_p1-1,vran+3*(n_ap-n_p1-1));
  }

  for (int i=1;i<n_ap;i++) p[p_perm[i]]=m_q[i];
  CE.Isotropic2Momenta(m_q[n_b],0.,0.,p[n_d1],p[n_d2],vran[rannum-2],vran[rannum-1]);
}

void VHAAG_res::CalculateS0(Cut_Data * cuts) 
{
  m_s0=0.; return;
}

int VHAAG_res::OType()
{
  return (1<<m_type);
}
