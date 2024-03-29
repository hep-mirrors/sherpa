#ifndef PHASIC_Channels_VHAAG_Threshold_h
#define PHASIC_Channels_VHAAG_Threshold_h

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <map>

namespace PHASIC {
  class VHAAG_Threshold : public Single_Channel {
    int      n_p1,m_type,n_d1,n_d2,n_b,n_ap;
    int      *p_perm;
    double   m_s0, m_thmass;
    ATOOLS::Vec4D* m_q;
    double*  p_s;
    Vegas* p_vegas;
    bool m_ownvegas, m_first;
    
    std::map<int,Vegas*>* p_sharedvegaslist;

    double PiFunc(double a1,double a2,
		  double s1b,double s2b,double c);
    void Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
	       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,int,double *ran);
    void Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,int,double *ran);
    void SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    void SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    void SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    double SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		       ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,int,double *ran);
    double Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,int,double *ran);
    double SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D& Q,
			     ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    double SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D& Q,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    double SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    void GenerateBosonMass(ATOOLS::Vec4D *p,double *ran);
    double BosonWeight(ATOOLS::Vec4D *p,double *ran);
    void GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			ATOOLS::Vec4D* q,double*,int n,double *ran);
    double BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
			ATOOLS::Vec4D* q,double*,int n,double *ran);
    void ConstructMomenta(double a1,double a2,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void CalculateS0(Cut_Data *);
     
    void Initialize(std::vector<int> perm, VHAAG_Threshold* ovl);

  public:

    VHAAG_Threshold(int _nin,int _nout,int pn,int d1,int d2,double th,VHAAG_Threshold* ovl);

    ~VHAAG_Threshold();

    void AddPoint(double Value);
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    void   MPISync()                 { p_vegas->MPISync(); }
    void   Optimize()                { p_vegas->Optimize(); }
    void   EndOptimize()             { p_vegas->EndOptimize(); }
    void   WriteOut(std::string pId) { if (m_ownvegas) p_vegas->WriteOut(pId); }
    void   ReadIn(std::string pId)   { if (m_ownvegas) p_vegas->ReadIn(pId); }

    int    RTH()                     { return m_thmass; }
    int    RD1()                     { return n_d1; }
    int    RD2()                     { return n_d2; }
    int    Type()                    { return m_type; }
    int    OType(); 
    std::map<int,Vegas*>* GetSharedVegasList() { return p_sharedvegaslist; }

    bool   OptimizationFinished()  { return p_vegas->Finished(); }
    void   ISRInfo(std::vector<int> &ts,std::vector<double> &ms,
		   std::vector<double> &ws) const;
  };
}
#endif

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
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

VHAAG_Threshold::VHAAG_Threshold(int _nin,int _nout,int pn,int d1,int d2,double th,VHAAG_Threshold* ovl)
{
  m_first=!ovl;
  m_nin=_nin; m_nout=_nout;n_ap=m_nin+m_nout-1;
  if (n_ap<5) {
    msg_Error()<<"Minimum number of final state particles for VHAAG_Threshold integrator is 4!"<<std::endl;
    Abort();
  }
  if (!ovl) {
    m_thmass=th;
    n_d1=d1;
    n_d2=d2;
    msg_Out()<<" Initialized HAAG with threshold m_{"<<n_d1<<","<<n_d2<<"} = "<<m_thmass<<" GeV."<<std::endl;
  }
  else {
    m_thmass = ovl->RTH();
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

void VHAAG_Threshold::Initialize(std::vector<int> perm, VHAAG_Threshold* ovl) 
{
  p_perm = new int[n_ap];
  p_perm[0] = 0;
  p_s = new double[n_ap];
  for (int i=0;i<n_ap;i++) p_s[i]=0.;
  msg_Tracking()<<"Init VHAAG_Threshold: 0";
  m_name = "VHAAG_Threshold";
  for (int i=1;i<n_ap;i++) {
    p_perm[i] = perm[i];
    if (perm[i]==1) n_p1=i;
    m_name+= "_";
    m_name+= std::to_string(p_perm[i]);
    msg_Tracking()<<" "<<p_perm[i];
  }

  m_rannum = 3*m_nout-4;;
  p_rans   = new double[m_rannum];
  m_q      = new Vec4D[n_ap];
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
      (*p_sharedvegaslist)[m_type] = new Vegas(m_rannum,100,Name());
//       (*p_sharedvegaslist)[m_type]->SetAutoOptimize(n_ap*10);//Min(m_nout,5)*15);
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
	(*p_sharedvegaslist)[m_type]->ConstChannel(m_rannum-4);
      } 
      m_ownvegas = true;
    } 
    p_vegas = (*p_sharedvegaslist)[m_type];
  } 
  else  {
    m_ownvegas = true;
    p_vegas = new Vegas(m_rannum,100,Name());
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
      p_vegas->ConstChannel(m_rannum-4);
    } 
  } 

  m_s0=-1.;
}

VHAAG_Threshold::~VHAAG_Threshold()
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

void VHAAG_Threshold::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,p_rans);
}

double VHAAG_Threshold::PiFunc(double a1,double a2,
		     double s1b,double s2b,double c)
{
  return 4.*(1.-c*c)*((1.-a2+s2b-s1b)*a2-s2b)
    -sqr((1.-2.*a1+s1b-s2b)+(1.-2.*a2-s1b+s2b)*c);
}

void VHAAG_Threshold::Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

void VHAAG_Threshold::Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

void VHAAG_Threshold::SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

void VHAAG_Threshold::SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

void VHAAG_Threshold::SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double s2,double *ran)
{
  double s = Q.Abs2();

  double a1min = 0.; 
  double a1max = 1.-s2/s; 
  double a1 = CE.AntennaMomenta(a1min,a1max,ran[0]);

  double phi=2.*M_PI*ran[1];

  ConstructMomenta(a1,phi,0.,s2,s,q1,p1,p2);
}

void VHAAG_Threshold::ConstructMomenta(double a1,double a2,
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

void VHAAG_Threshold::ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)
{
  double ps = 0.25*(sqr(s-s1-s2)-4.*s1*s2)/s;
  if (q1.PPerp()!=0.||!IsZero(q1.Abs2()/Max(1.,sqr(q1[0])),1.0e-6)) {
    msg_Error()<<" Error in"<<std::endl
	       <<"ConstructMomenta(double a1,double phi,double s1,double s2,double s,"<<std::endl
	       <<"                 ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)!"<<std::endl
	       <<" q1 must be in beam direction and massless!   q1="<<q1<<" ("<<q1.Abs2()<<")"<<std::endl;
    Abort();
  }
  Vec3D e1  = Vec3D(q1)/q1[0];
  double v1 = sqrt(ps+s1)-a1*sqrt(s);

  double cc = sqrt(ps-v1*v1); 
  Vec3D pv(cc*cos(phi),cc*sin(phi),v1*e1[3]);

  p1 = Vec4D(sqrt(ps+s1),pv);
  p2 = Vec4D(sqrt(ps+s2),-1.*pv);
}
    
double VHAAG_Threshold::SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

double VHAAG_Threshold::Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

double VHAAG_Threshold::SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
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

double VHAAG_Threshold::SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
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

double VHAAG_Threshold::SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

void VHAAG_Threshold::GenerateBosonMass(ATOOLS::Vec4D *p,double *ran)
{
  double smax=(p[0]+p[1]).Abs2();
  p_s[n_b]=CE.ThresholdMomenta(1.0,m_thmass,0.,smax,ran[m_rannum-3]);
}

double VHAAG_Threshold::BosonWeight(ATOOLS::Vec4D *p,double *ran)
{
  double smax=(p[0]+p[1]).Abs2();
  double w=CE.ThresholdWeight(1.0,m_thmass,0.,smax,m_q[n_b].Abs2(),ran[m_rannum-3]);
  w*=CE.Isotropic2Weight(p[n_d1],p[n_d2],ran[m_rannum-2],ran[m_rannum-1]);
  return w;
}


void VHAAG_Threshold::GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
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

double VHAAG_Threshold::BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
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

void VHAAG_Threshold::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  CalculateS0(cuts);
  double wt=1.;

  if (n_ap==4) {
    Vec4D Q(p[0]+p[1]);
    wt=SingleSplitF0Weight(p[0],Q,p[2],p[3],p_s[3],p_rans);  
    m_weight = p_vegas->GenerateWeight(p_rans)/wt/pow(2.*M_PI,2);
    return;
  }

  for (int i=0;i<n_ap;i++) m_q[i]=p[p_perm[i]];
  m_q[n_b]=p[n_d1]+p[n_d2];
  p_s[n_b]=m_q[n_b].Abs2();
  Vec4D Q(m_q[0]+m_q[n_p1]);
  
  if (n_p1==1){
    Vec4D P;
    wt*=BranchWeight(m_q[2],P,&(m_q[3]),&(p_s[3]),n_ap-3,p_rans+3);
    wt*=Split0Weight(m_q[1],Q,m_q[2],P,2,3,p_rans);    
  }
  else if (n_p1==n_ap-1){
    Vec4D P;
    wt*=BranchWeight(m_q[1],P,&(m_q[2]),&(p_s[2]),n_ap-3,p_rans+3);
    wt*=Split0Weight(m_q[0],Q,m_q[1],P,1,2,p_rans);    
  }
  else if (n_p1==2){
    Vec4D P;
    wt*=BranchWeight(m_q[2],P,&(m_q[3]),&(p_s[3]),n_ap-3,p_rans+3);
    wt*=Split0Weight(m_q[0],Q,m_q[1],P,1,3,p_rans);    
  }
  else if (n_p1==n_ap-2){
    Vec4D P;
    wt*=BranchWeight(m_q[0],P,&(m_q[1]),&(p_s[1]),n_ap-3,p_rans+3);
    wt*=Split0Weight(m_q[n_p1],Q,m_q[n_ap-1],P,n_ap-1,1,p_rans);    
  }
  else if (n_p1<=(n_ap-1)/2) {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[0],Q1,&(m_q[1]),&(p_s[1]),n_p1-1,p_rans+4);
    wt*=BranchWeight(m_q[n_p1],Q2,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,p_rans+3*(n_p1-1));
    wt*=SplitWeight(m_q[0],Q,Q1,Q2,0,n_p1,p_rans);
  }
  else {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[n_p1],Q1,&(m_q[n_p1+1]),&(p_s[n_p1+1]),n_ap-n_p1-1,p_rans+4);
    wt*=BranchWeight(m_q[0],Q2,&(m_q[1]),&(p_s[1]),n_p1-1,p_rans+3*(n_ap-n_p1-1));
    wt*=SplitWeight(m_q[n_p1],Q,Q1,Q2,n_p1,0,p_rans);
  }
  wt*=BosonWeight(p,p_rans);
  double vw = p_vegas->GenerateWeight(p_rans);
  m_weight = vw/wt/pow(2.*M_PI,m_nout*3.-4.);
}

void VHAAG_Threshold::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *ran)
{
  CalculateS0(cuts);
  double *vran = p_vegas->GeneratePoint(ran);
  for(int i=0;i<m_rannum;i++) p_rans[i]=vran[i];

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
  CE.Isotropic2Momenta(m_q[n_b],0.,0.,p[n_d1],p[n_d2],vran[m_rannum-2],vran[m_rannum-1]);
}

void VHAAG_Threshold::CalculateS0(Cut_Data * cuts) 
{
  m_s0=0.; return;
}

int VHAAG_Threshold::OType()
{
  return (1<<m_type);
}

void VHAAG_Threshold::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,std::vector<double> &ws) const
{
  if (!m_first) return;
  ts.push_back(2);
  ms.push_back(m_thmass);
  ws.push_back(0.0);
}

namespace PHASIC {

  class VHAAG_Threshold_Channel_Generator: public Channel_Generator {
  private:

    int m_d1, m_d2;
    double m_th;

  public:
    
    VHAAG_Threshold_Channel_Generator(const Channel_Generator_Key &key):
      Channel_Generator(key), m_d1(2), m_d2(3), m_th(40.0)
    {
      size_t bpos(key.m_key.find('[')), epos(key.m_key.rfind(']'));
      if (bpos==std::string::npos || epos==std::string::npos) {
        Scoped_Settings s{ Settings::GetMainSettings()["VHAAG"] };
	m_th=s["TH"].SetDefault(40.0).Get<int>();
	m_d1=s["D1"].SetDefault(2).Get<int>();
	m_d2=s["D2"].SetDefault(3).Get<int>();
	return;
      }
      Data_Reader read(":",",","#","=");
      read.SetString(key.m_key.substr(bpos+1,epos-bpos-1));
      if (!read.ReadFromString(m_d1,"I")) m_d1=2;
      if (!read.ReadFromString(m_d2,"J")) m_d2=3;
      if (!read.ReadFromString(m_th,"T")) m_th=40.0;
    }

    int GenerateChannels()
    {
      int m_nin=p_proc->NIn(), m_nout=p_proc->NOut();
      VHAAG_Threshold *firsthaag=NULL,*hlp=NULL;
      Permutation pp(m_nin+m_nout-3);
      for (int j=0;j<pp.MaxNumber();j++) {
	p_mc->Add(hlp=new VHAAG_Threshold(m_nin,m_nout,2*j,m_d1,m_d2,m_th,firsthaag));
	if (!firsthaag) firsthaag=hlp;
	p_mc->Add(hlp=new VHAAG_Threshold(m_nin,m_nout,2*j+1,m_d1,m_d2,m_th,firsthaag));
	if (!firsthaag) firsthaag=hlp;
      }
      return 0;
    }

  };// end of class VHAAG_Threshold_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(VHAAG_Threshold_Channel_Generator,"VHAAG_Threshold",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,VHAAG_Threshold_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new VHAAG_Threshold_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    VHAAG_Threshold_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Vegas-improved HAAG integrator for resonance + jets";
}
