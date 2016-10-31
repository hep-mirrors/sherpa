#include "PHASIC++/Channels/VHAAG_ND.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"
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
const double s_xx  =.3;
const bool ap =false;
const bool apF=false;

VHAAG_ND::VHAAG_ND(int _nin,int _nout,int pn, VHAAG_ND* ovl) 
{
  Permutation pp(_nin+_nout-1);
  int* tp=pp.Get(pn);
  std::vector<int> perm(_nin+_nout);
  perm[0] = 0;
  for (int i=1;i<_nin+_nout;i++) perm[i] = tp[i-1]+1;
  Initialize(_nin,_nout,perm,ovl);
}

void VHAAG_ND::Initialize(int _nin,int _nout,std::vector<int> perm, VHAAG_ND* ovl) 
{
  nin=_nin; nout=_nout;
  Permutation pp(nin+nout-1);
  p_perm = new int[nin+nout];
  p_mrep = new int[nin+nout];
  p_perm[0] = p_mrep[0] = 0;
  msg_Tracking()<<"Init VHAAG_ND: 0";
  name = "VHAAG_ND";
  char hlp[4];
  for (int i=1;i<nin+nout;i++) {
    p_perm[i] = perm[i];
    p_mrep[p_perm[i]] = i;
    if (perm[i]==1) n_p1=i;
    name+= "_";
    sprintf(hlp,"%i",p_perm[i]);
    name+= std::string(hlp);
    msg_Tracking()<<" "<<p_perm[i];
  }

  if (n_p1==2||n_p1==nout) {
    msg_Error()<<"Invalid channel for VHAAG_ND!"<<std::endl;
    Abort();
  }
// std::cout<<" n_p1="<<n_p1<<std::endl;

  rannum = 3*nout-4;;
  rans  = new double[rannum];
  m_q = new Vec4D[nin+nout];
  m_ownvegas = false;
  if (ovl) p_sharedvegaslist = ovl->GetSharedVegasList();
  else p_sharedvegaslist = NULL;
  if (p_sharedvegaslist==NULL) {
    p_sharedvegaslist = new Vegas*[nout];
    for (int i=0;i<nout;i++) p_sharedvegaslist[i]=NULL;
  }

  m_type=Min(n_p1-1,nout-(n_p1-1));
  msg_Tracking()<<" n_p1="<<n_p1<<" type="<<m_type<<std::endl;
  int vs=m_type;

  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(rpa->GetPath());
  dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  int size=dr.GetValue<std::string>("VHAAG_AUTOOPT","On")=="On"?1:2;
#ifdef USING__MPI
  size=MPI::COMM_WORLD.Get_size();
#endif
  if (1) {
    if (p_sharedvegaslist[vs]==NULL) {
      p_sharedvegaslist[vs] = new Vegas(rannum,100,Name());
      if (size<2) p_sharedvegaslist[vs]->SetAutoOptimize(Min(nout,5)*100);
      if (1) {
	if (m_type<2) {
	  for (int i=0;i<=nout-3;i++) p_sharedvegaslist[vs]->ConstChannel(2+3*i);
	} 
	else {
	  for (int i=0;i<=m_type-2;i++) p_sharedvegaslist[vs]->ConstChannel(3+3*i);
	  p_sharedvegaslist[vs]->ConstChannel(3*(m_type-1)+2);
	  for (int i=0;i<nout-m_type-2;i++) p_sharedvegaslist[vs]->ConstChannel(2+3*m_type+3*i);
	} 
	p_sharedvegaslist[vs]->ConstChannel(rannum-1);
      } 
      m_ownvegas = true;
    } 
    p_vegas = p_sharedvegaslist[vs];
  } 
  else  {
    m_ownvegas = true;
    p_vegas = new Vegas(rannum,100,Name());
    if (size<2) p_vegas->SetAutoOptimize(500);
  } 

  m_s0=-1.;
}

VHAAG_ND::~VHAAG_ND()
{
  delete[] p_perm;
  delete[] p_mrep;
  delete[] m_q;
  if (m_ownvegas) {
    delete p_vegas;
    if (p_sharedvegaslist) p_sharedvegaslist[m_type]=0;
  }
  if (p_sharedvegaslist) {
    int empty=1;
    for (int i=0;i<nout;i++) if (p_sharedvegaslist[i]!=0) empty=0;
    if (empty) delete[] p_sharedvegaslist;
  }
}

void VHAAG_ND::MPISync()
{
#ifdef USING__MPI
  if (m_ownvegas) p_vegas->MPISync();
#endif
}

void VHAAG_ND::Optimize()
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_size()<2) return;
  if (m_ownvegas) p_vegas->Optimize();
#endif
}

void VHAAG_ND::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}

double VHAAG_ND::PiFunc(double a1,double a2,
		     double s1b,double s2b,double c)
{
  return 4.*(1.-c*c)*((1.-a2+s2b-s1b)*a2-s2b)
    -sqr((1.-2.*a1+s1b-s2b)+(1.-2.*a2-s1b+s2b)*c);
}

void VHAAG_ND::Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
	       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n1, int n2,double *ran)
{
  double s = (q1+q2).Abs2();
  double s1min = double((n1-1)*(n1-2)/2)*m_s0;
  double s2min = double((n2-1)*(n2-2)/2)*m_s0;
  double s1max = sqr(sqrt(s)-sqrt(s2min));
  double exp = s_sp1; if (s1min==0.) exp = s_xx;
  double s1 = CE.MasslessPropMomenta(exp,s1min,s1max,ran[0]);

  exp = s_sp2; if (s2min==0.) exp = s_xx;
  double s2max = sqr(sqrt(s)-sqrt(s1));
  double s2 = CE.MasslessPropMomenta(exp,s2min,s2max,ran[1]);

  double pb0 =0.5*(s+s1-s2)/s;
  double pb  =sqrt(pb0*pb0-s1/s);

  double a1min = pb0-pb;//Max(0.5*m_s0*(n1-1)/(q1*q2),pb0-pb); 
  double a1max = pb0+pb;//Min(1.-0.5*m_s0*(n2-1)/(q1*q2),pb0+pb); //to check!!!
  exp = a1_sp; if (a1min==0.) exp = s_xx;
  double a1 = CE.MasslessPropMomenta(exp,a1min,a1max,ran[2]);

  double phi=2.*M_PI*ran[3];

  ConstructMomenta(a1,phi,s1,s2,s,q1,p1,p2);
//     std::cout<<"Split generated: "<<p1<<p2<<std::endl;
//     std::cout<<" s1: "<<s1min<<" < "<<s1<<" < "<<s1max<<std::endl;
//     std::cout<<" s2: "<<s2min<<" < "<<s2<<" < "<<s2max<<std::endl;
//     std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<" "<<a1*(q1*q2)<<std::endl;
}

void VHAAG_ND::Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n2,double *ran)
{
  double s = (q1+q2).Abs2();
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double exp = s_s; if (smin==0.) exp = s_xx;
  double s2 = CE.MasslessPropMomenta(exp,smin,smax,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = 0.;//Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = 2.*p;//Min(1.-a1min*n2,2.*p); 
  double a1 = CE.MasslessPropMomenta(s_xx,a1min,a1max,ran[1]);

  double phi=2.*M_PI*ran[2];

  ConstructMomenta(a1,phi,0.,s2,s,q1,p1,p2);
//    std::cout<<"Split0 generated: "<<p1<<p2<<std::endl;
//    std::cout<<" s2: "<<smin<<" < "<<s2<<" < "<<smax<<std::endl;
//    std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<" "<<a1*(q1*q2)<<std::endl;
}

void VHAAG_ND::SplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n2,double *ran)
{
  double s = (q1+q2).Abs2();
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double exp = 1.; if (smin==0.) exp = s_xx;
  double s2 = CE.MasslessPropMomenta(exp,smin,smax,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = 0.;//Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = Min(1.-a1min*n2,2.*p); 
  double a1 = a1min+(a1max-a1min)*ran[1];

  double phi=2.*M_PI*ran[2];

  ConstructMomenta(a1,phi,0.,s2,s,q1,p1,p2);
//    std::cout<<"Split0 generated: "<<p1<<p2<<std::endl;
//    std::cout<<" s2: "<<smin<<" < "<<s2<<" < "<<smax<<std::endl;
//    std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<" "<<a1*(q1*q2)<<std::endl;
}

void VHAAG_ND::Split1(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n2,double *ran)
{
  double s = (q1+q2).Abs2();
  double smin = double(n2*(n2-1)/2)*m_s0;
  double smax = Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double s2 = CE.MasslessPropMomenta(s_s,smin,smax,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = Min(1.-a1min*n2,2.*p);
  double a1 = CE.AntennaMomenta(a1min,a1max,ran[1]);

  double phi=2.*M_PI*ran[2];

  ConstructMomenta(a1,phi,0.,s2,s,q1,p1,p2);
//   std::cout<<"Split1 generated: "<<p1<<p2<<std::endl;
//   std::cout<<" s2: "<<smin<<" < "<<s2<<" < "<<smax<<std::endl;
//   std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<" "<<a1*(q1*q2)<<std::endl;
}

void VHAAG_ND::SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
		     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int n2,double *ran)
{
  Poincare qb(Q);
  qb.Boost(q1);
  qb.Boost(q2);
  double s = Q.Abs2();
  double v  = Vec3D(q1)*Vec3D(q2)/(q1[0]*q2[0]);
  double a1min = 0.;//0.5*m_s0/(sqrt(s)*q1[0]);
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s*(1.-a1min));
  double exp = s_i; if (smin==0.) exp = s_xx;
  double s2 = CE.MasslessPropMomenta(exp,smin,smax,ran[0]);

  double a1max = 1.-s2/s;//Min(1.-a1min*n2,1.-s2/s); 
  exp = a1_i; if (a1min==0.) exp = s_xx;
  double a1 = CE.MasslessPropMomenta(exp,a1min,a1max,ran[1]);

  if (ap) {
    double hlp0 = 0.5*(1.+s2/s+v*(1.-s2/s-2.*a1));
    double hlp = sqrt((1.-v*v)*(1.-s2/s-a1)*a1);
    double a2max = Min(hlp0+hlp,1.-0.5*m_s0/(sqrt(s)*q2[0]));
    double a2min = Max(hlp0-hlp,0.5*m_s0*n2/(sqrt(s)*q2[0]));
    if (a2max<=a2min) {
      a2min=Max(0.,hlp0-hlp);  
      a2max=Min(1.,hlp0+hlp);  
    }
    double a2 = CE.MasslessPropMomenta(a2_i,a2min,a2max,ran[2]);
    
    ConstructMomenta(a1,a2,0.,s2,s,q1,q2,p1,p2);
  }
  else {
    double phi=2.*M_PI*ran[2];
    Vec4D qq(1.,0.,0.,1.);
    ConstructMomenta(a1,phi,0.,s2,s,qq,p1,p2);
    Poincare rot(qq,q1);
    rot.Rotate(p1);
    rot.Rotate(p2);
  }
  qb.BoostBack(p1);
  qb.BoostBack(p2);
//   std::cout<<"SingleSplit generated: "<<p1<<p2<<std::endl;
//   std::cout<<q1<<q2<<Q<<v<<std::endl;
//   std::cout<<" s2: "<<smin<<" < "<<s2<<" < "<<smax<<std::endl;
//   std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<std::endl;
//   std::cout<<" a2: "<<a2min<<" < "<<a2<<" < "<<a2max<<std::endl;
}

void VHAAG_ND::SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
		      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double *ran)
{
  Poincare qb(Q);
  qb.Boost(q1);
  qb.Boost(q2);
  double s = Q.Abs2();
  double v  = Vec3D(q1)*Vec3D(q2)/(q1[0]*q2[0]);

  double a1min = 0.;//0.5*m_s0/(sqrt(s)*q1[0]); 
  double a1max = 1.-a1min; 
  double a1 = CE.MasslessPropMomenta(s_xx,a1min,a1max,ran[0]);

  if (apF) {
    double hlp0 = 0.5*(1.+v*(1.-2.*a1));
    double hlp = sqrt((1.-v*v)*(1.-a1)*a1);
    double a2max = Min(hlp0+hlp,1.-0.5*m_s0/(sqrt(s)*q2[0]));
    double a2min = Max(hlp0-hlp,0.5*m_s0/(sqrt(s)*q2[0]));
    if (a2max<=a2min) {
      a2min=Max(0.,hlp0-hlp);  
      a2max=Min(1.,hlp0+hlp);  
    }
    double a2 = CE.MasslessPropMomenta(a2_iF,a2min,a2max,ran[1]);
    
    ConstructMomenta(a1,a2,0.,0.,s,q1,q2,p1,p2);
  }
  else {
    double phi=2.*M_PI*ran[1];
    Vec4D qq(1.,0.,0.,1.);
    ConstructMomenta(a1,phi,0.,0.,s,qq,p1,p2);
    Poincare rot(qq,q1);
    rot.Rotate(p1);
    rot.Rotate(p2);
  }

  qb.BoostBack(p1);
  qb.BoostBack(p2);
//     std::cout<<"SingleSplitF generated: "<<p1<<p2<<std::endl;
//     std::cout<<" a1: "<<a1min<<" < "<<a1<<" < "<<a1max<<std::endl;
//     std::cout<<" a2: "<<a2min<<" < "<<a2<<" < "<<a2max<<std::endl;
}

void VHAAG_ND::SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double *ran)
{
  double s = (q1+q2).Abs2();

  double a1min = 0.5*m_s0/(q1*q2); 
  double a1max = 1.-a1min; 
  double a1 = CE.AntennaMomenta(a1min,a1max,ran[0]);

  double phi=2.*M_PI*ran[1];

  ConstructMomenta(a1,phi,0.,0.,s,q1,p1,p2);
}

void VHAAG_ND::ConstructMomenta(double a1,double a2,
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

void VHAAG_ND::ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)
{
  double ps = 0.25*(sqr(s-s1-s2)-4.*s1*s2)/s;
  if (q1.PPerp()!=0.) {
    msg_Error()<<" Error in"<<std::endl
	       <<"ConstructMomenta(double a1,double phi,double s1,double s2,double s,"<<std::endl
	       <<"                 ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2)!"<<std::endl
	       <<" q1 must be in beam direction!   q1="<<q1<<std::endl;
    Abort();
  }
  Vec3D e1  = Vec3D(q1)/q1[0];
  double v1 = sqrt(ps+s1)-a1*sqrt(s);

  double cc = sqrt(ps-v1*v1); 
  Vec3D pv(cc*cos(phi),cc*sin(phi),v1*e1[3]);

  p1 = Vec4D(sqrt(ps+s1),pv);
  p2 = Vec4D(sqrt(ps+s2),-1.*pv);
}
    
double VHAAG_ND::SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n1, int n2,double *ran)
{
  double wt=1.;
  double s = (q1+q2).Abs2();
  double s1min = double((n1-1)*(n1-2)/2)*m_s0;
  double s2min = double((n2-1)*(n2-2)/2)*m_s0;
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

  s1 = 1./(q1*q2);
  double a1min = pb0-pb;//Max(0.5*m_s0*(n1-1)/(q1*q2),pb0-pb); 
  double a1max = pb0+pb;//Min(1.-0.5*m_s0*(n2-1)/(q1*q2),pb0+pb); //to check!!!
  exp = a1_sp; if (a1min==0.) exp = s_xx;
  double a1 = (q1*p1)*s1;
  wt*= CE.MasslessPropWeight(exp,a1min,a1max,a1,ran[2]);

  wt*= 2./M_PI;
  ran[3]=p1.Phi()/(2.*M_PI);
  if(ran[3]<0.) ran[3]+=1.;
  return wt;
}

double VHAAG_ND::SplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			   ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n2,double *ran)
{
  double wt=1.;
  double s = (q1+q2).Abs2();
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double exp = 1.; if (smin==0.) exp = s_xx;
  double s2 = p2.Abs2(); 
  wt*= CE.MasslessPropWeight(exp,smin,smax,s2,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = Min(1.-a1min*n2,2.*p); 
  double a1 = q1*p1/(q1*q2);
  wt*= 1./(a1max-a1min);
  ran[1] = (a1-a1min)/(a1max-a1min);

  wt*= 2./M_PI;
  ran[2]=p1.Phi()/(2.*M_PI);
  if(ran[2]<0.) ran[2]+=1.;
  return wt;
}

double VHAAG_ND::Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			   ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n2,double *ran)
{
  double wt=1.;
  double s = (q1+q2).Abs2();
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double exp = s_s; if (smin==0.) exp = s_xx;
  double s2 = p2.Abs2(); 
  wt*= CE.MasslessPropWeight(exp,smin,smax,s2,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = 0.;//Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = 2.*p;//Min(1.-a1min*n2,2.*p); 
  double a1 = q1*p1/(q1*q2);
  wt*= CE.MasslessPropWeight(s_xx,a1min,a1max,a1,ran[1]);

  wt*= 2./M_PI;
  ran[2]=p1.Phi()/(2.*M_PI);
  if(ran[2]<0.) ran[2]+=1.;
  return wt;
}

double VHAAG_ND::Split1Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			 ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n2,double *ran)
{
  double wt=1.;
  double s = (q1+q2).Abs2();
  double smin = double(n2*(n2-1)/2)*m_s0;
  double smax = Min(s-n2*m_s0,s-2.*sqrt(s*m_s0));
  double s2 = p2.Abs2(); 
  wt*= CE.MasslessPropWeight(s_s,smin,smax,s2,ran[0]);

  double p = 0.5*(s-s2)/s;
  double a1min = Max(0.5*m_s0/(q1*q2),p*(1.-sqrt(1.-m_s0/(p*p*s)))); 
  double a1max = Min(1.-a1min*n2,2.*p); 
  double a1 = q1*p1/(q1*q2);
  wt*= CE.AntennaWeight(a1min,a1max,a1,ran[1]);

  wt*= 2./M_PI;
  ran[2]=p1.Phi()/(2.*M_PI);
  if(ran[2]<0.) ran[2]+=1.;
  return wt;
}

double VHAAG_ND::SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D &Q,
				ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int n2,double *ran)
{
  double wt=1.;
  Q = p1+p2;
  double s = Q.Abs2();
  double s2 = p2.Abs2();
  double a1min = 0.;//0.5*m_s0/(sqrt(s)*q1[0]);
  double smin = double((n2-1)*(n2-2)/2)*m_s0;
  double smax = s;//Min(s-n2*m_s0,s*(1.-a1min));
  double exp = s_i; if (smin==0.) exp = s_xx;
  wt*= CE.MasslessPropWeight(exp,smin,smax,s2,ran[0]);

  double a1max = 1.-s2/s;//Min(1.-a1min*n2,1.-s2/s); 
  exp = a1_i; if (a1min==0.) exp = s_xx;
  double a1 = q1*p1/(q1*Q);
  double a2 = q2*p2/(q2*Q);
  wt*= CE.MasslessPropWeight(exp,a1min,a1max,a1,ran[1]);

  Poincare qb(Q);
  qb.Boost(q1);
  qb.Boost(q2);
  if (ap) {
    double v  = Vec3D(q1)*Vec3D(q2)/(q1[0]*q2[0]);
    double hlp0 = 0.5*(1.+s2/s+v*(1.-s2/s-2.*a1));
    double hlp = sqrt((1.-v*v)*(1.-s2/s-a1)*a1);
    double a2max = Min(hlp0+hlp,1.-0.5*m_s0/(sqrt(s)*q2[0]));
    double a2min = Max(hlp0-hlp,0.5*m_s0*n2/(sqrt(s)*q2[0]));
    if (a2max<=a2min) {
      a2min=Max(0.,hlp0-hlp);  
      a2max=Min(1.,hlp0+hlp);  
    }
    wt*= CE.MasslessPropWeight(a2_i,a2min,a2max,a2,ran[2]);
    
    hlp= PiFunc(a1,a2,0.,s2/s,v);
    if (hlp>0.) wt*=sqrt(hlp);
    else {
      msg_Error()<<"Error in VHAAG_ND::SingleSplitFWeight!"<<std::endl
		 <<"PiFunc("<<a1<<","<<a2<<",0,0,"<<v<<") = "<<hlp<<std::endl;
      wt=0.;
    }
  }
  else {
    wt*=2./M_PI;
    Vec4D qq(1.,0.,0.,1.);
    qb.Boost(p1);
    Poincare rot(qq,q1);
    rot.RotateBack(p1);
    ran[2]=p1.Phi()/(2.*M_PI);
    if(ran[2]<0.) ran[2]+=1.;
  }
  return wt;
}

double VHAAG_ND::SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D &Q,
				 ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double *ran)
{
  double wt=1.;
  Q = p1+p2;
  double s=Q.Abs2();
  double a1min = 0.;//0.5*m_s0/(sqrt(s)*q1[0]); 
  double a1max = 1.-a1min; 
  double a1 = q1*p1/(q1*Q);
  double a2 = q2*p2/(q2*Q);
  wt*= CE.MasslessPropWeight(s_xx,a1min,a1max,a1,ran[0]);

  Poincare qb(Q);
  qb.Boost(q1);
  qb.Boost(q2);
  if (apF) {
    double v  = Vec3D(q1)*Vec3D(q2)/(q1[0]*q2[0]);
    double hlp0 = 0.5*(1.+v*(1.-2.*a1));
    double hlp = sqrt((1.-v*v)*(1.-a1)*a1);
    double a2x,a2n;
    double a2max = a2x=Min(hlp0+hlp,1.-0.5*m_s0/(sqrt(s)*q2[0]));
    double a2min = a2n=Max(hlp0-hlp,0.5*m_s0/(sqrt(s)*q2[0]));
    if (a2max<=a2min) {
      a2min=Max(0.,hlp0-hlp);  
      a2max=Min(1.,hlp0+hlp);  
    }
    wt*= CE.MasslessPropWeight(a2_iF,a2min,a2max,a2,ran[1]);
    
    hlp= PiFunc(a1,a2,0.,0.,v);
    if (hlp>0.) wt*=sqrt(hlp);
    else {
      msg_Error()<<"Error in VHAAG_ND::SingleSplitFWeight!"<<std::endl
		 <<"PiFunc("<<a1<<","<<a2<<",0,0,"<<v<<") = "<<hlp<<std::endl;
      wt=-1.;
    }
  }
  else {
    wt*=2./M_PI;
    Vec4D qq(1.,0.,0.,1.);
    qb.Boost(p1);
    Poincare rot(qq,q1);
    rot.RotateBack(p1);
    ran[1]=p1.Phi()/(2.*M_PI);
    if(ran[1]<0.) ran[1]+=1.;
  }
  return wt;
}

double VHAAG_ND::SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
				 ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double *ran)
{
  double wt=1.;
  double a1min = 0.5*m_s0/(q1*q2); 
  double a1max = 1.-a1min; 
  double a1 = q1*p1/(q1*q2);
  wt*= CE.AntennaWeight(a1min,a1max,a1,ran[0]);

  wt*=2./M_PI;
  ran[1]=p1.Phi()/(2.*M_PI);
  if(ran[1]<0.) ran[1]+=1.;

  return wt;
}

void VHAAG_ND::GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
			   ATOOLS::Vec4D* q,int n,double *ran)
{
  Vec4D r1 = q1;
  Vec4D rQ = Q;
  for (int i=n;i>2;i--) {
    SingleSplit(r1,q2,rQ,q[n-i],rQ,i-1,ran);
    ran+=3;
    r1 = q[n-i];
  }
  SingleSplitF(r1,q2,rQ,q[n-2],q[n-1],ran);
}

double VHAAG_ND::BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D &Q,
			 ATOOLS::Vec4D* q,int n,double *ran)
{
  double wt=1.;
  ran+= 3*(n-2);
  wt*=SingleSplitFWeight(q[n-3],q2,Q,q[n-2],q[n-1],ran);

  for (int i=3;i<=n;i++) {
    ran-= 3;
    wt*=SingleSplitWeight(q[n-i-1],q2,Q,q[n-i],Q,i-1,ran);
  }
  return wt;
}

void VHAAG_ND::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  CalculateS0(cuts);
  double wt=1.;

  if (nout==2) {
    wt=SingleSplitF0Weight(p[0],p[1],p[2],p[3],rans);  
    weight = p_vegas->GenerateWeight(rans)/wt/pow(2.*M_PI,2);
    return;
  }

  for (int i=0;i<nin+nout;i++) m_q[i]=p[p_perm[i]];
  
  if (n_p1==1){
    Vec4D Q;
    wt*=BranchWeight(m_q[2],m_q[0],Q,&(m_q[3]),nout-1,rans+3);
    wt*=Split0Weight(m_q[1],m_q[0],m_q[2],Q,nout-1,rans);    
  }
  else if (n_p1==nout+1){
    Vec4D Q;
    wt*=BranchWeight(m_q[1],m_q[n_p1],Q,&(m_q[2]),nout-1,rans+3);
    wt*=Split0Weight(m_q[0],m_q[n_p1],m_q[1],Q,nout-1,rans);    
  }
  else if (n_p1==2){
    Vec4D Q;
    wt*=BranchWeight(m_q[2],m_q[0],Q,&(m_q[3]),nout-1,rans+3);
    wt*=Split1Weight(m_q[0],m_q[2],m_q[1],Q,nout-1,rans);    
  }
  else if (n_p1==nout){
    Vec4D Q;
    wt*=BranchWeight(m_q[0],m_q[n_p1],Q,&(m_q[1]),nout-1,rans+3);
    wt*=Split1Weight(m_q[n_p1],m_q[0],m_q[nout+1],Q,nout-1,rans);    
  }
  else if (n_p1<=(nout+1)/2) {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[0],m_q[n_p1],Q1,&(m_q[1]),n_p1-1,rans+4);
    wt*=BranchWeight(m_q[n_p1],m_q[0],Q2,&(m_q[n_p1+1]),nout-n_p1+1,rans+3*(n_p1-1));
    wt*=SplitWeight(m_q[0],m_q[n_p1],Q1,Q2,n_p1-1,nout-n_p1+1,rans);
  }
  else {
    Vec4D Q1,Q2;
    wt*=BranchWeight(m_q[n_p1],m_q[0],Q1,&(m_q[n_p1+1]),nout-n_p1+1,rans+4);
    wt*=BranchWeight(m_q[0],m_q[n_p1],Q2,&(m_q[1]),n_p1-1,rans+3*(nout-n_p1+1));
    wt*=SplitWeight(m_q[n_p1],m_q[0],Q1,Q2,nout-n_p1+1,n_p1-1,rans);
  }
  double vw = p_vegas->GenerateWeight(rans);
  weight = vw/wt/pow(2.*M_PI,nout*3.-4.);
}

void VHAAG_ND::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *ran)
{
  CalculateS0(cuts);
  double *vran = p_vegas->GeneratePoint(ran);
  for(int i=0;i<rannum;i++) rans[i]=vran[i];

  if (nout==2) {
    SingleSplitF0(p[0],p[1],p[2],p[3],ran);  
    return;
  }

  m_q[0] = p[0];
  m_q[n_p1] = p[1];

  if (n_p1==1){
    Vec4D Q;
    Split0(m_q[1],m_q[0],m_q[2],Q,nout-1,vran);    
    GenerateBranch(m_q[2],m_q[0],Q,&(m_q[3]),nout-1,vran+3);
  }
  else if (n_p1==nout+1){
    Vec4D Q;
    Split0(m_q[0],m_q[n_p1],m_q[1],Q,nout-1,vran);    
    GenerateBranch(m_q[1],m_q[n_p1],Q,&(m_q[2]),nout-1,vran+3);
  }
  else if (n_p1==2){
    Vec4D Q;
    Split1(m_q[0],m_q[2],m_q[1],Q,nout-1,vran);    
    GenerateBranch(m_q[2],m_q[0],Q,&(m_q[3]),nout-1,vran+3);
  }
  else if (n_p1==nout){
    Vec4D Q;
    Split1(m_q[n_p1],m_q[0],m_q[nout+1],Q,nout-1,vran);    
    GenerateBranch(m_q[0],m_q[n_p1],Q,&(m_q[1]),nout-1,vran+3);
  }
  else if (n_p1<=(nout+1)/2) {
    Vec4D Q1,Q2;
    Split(m_q[0],m_q[n_p1],Q1,Q2,n_p1-1,nout-n_p1+1,vran);
    GenerateBranch(m_q[0],m_q[n_p1],Q1,&(m_q[1]),n_p1-1,vran+4);
    GenerateBranch(m_q[n_p1],m_q[0],Q2,&(m_q[n_p1+1]),nout-n_p1+1,vran+3*(n_p1-1));
  }
  else {
    Vec4D Q1,Q2;
    Split(m_q[n_p1],m_q[0],Q1,Q2,nout-n_p1+1,n_p1-1,vran);
    GenerateBranch(m_q[n_p1],m_q[0],Q1,&(m_q[n_p1+1]),nout-n_p1+1,vran+4);
    GenerateBranch(m_q[0],m_q[n_p1],Q2,&(m_q[1]),n_p1-1,vran+3*(nout-n_p1+1));
  }

   for (int i=1;i<nin+nout;i++) p[p_perm[i]]=m_q[i];
}

void VHAAG_ND::CalculateS0(Cut_Data * cuts) 
{
  if (m_s0>0.) return;
  m_s0 = 0.;
  for (int i=0;i<cuts->ncut;i++) {
    for (int j=i+1;j<cuts->ncut;j++) {
      if (m_s0<cuts->scut[i][j]) m_s0 = cuts->scut[i][j];
    }
  }
}

int VHAAG_ND::OType()
{
  return (1<<m_type);
}

namespace PHASIC {

  class VHAAG_ND_Channel_Generator: public Channel_Generator {
  public:
    
    VHAAG_ND_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      int m_nin=p_proc->NIn(), m_nout=p_proc->NOut();
      if (m_nin==2 && m_nout==2) {
	p_mc->Add(new S1Channel(m_nin,m_nout,(Flavour*)&p_proc->Flavours().front()));
	p_mc->Add(new T1Channel(m_nin,m_nout,(Flavour*)&p_proc->Flavours().front()));
	p_mc->Add(new U1Channel(m_nin,m_nout,(Flavour*)&p_proc->Flavours().front()));
      }
      else {
	VHAAG_ND *firsthaag=NULL,*hlp=NULL;
	Permutation pp(m_nin+m_nout-1);
	for (int j=0;j<pp.MaxNumber();j++) {
	  int* pm = pp.Get(j);
	  if (pm[1]!=0&&pm[m_nin+m_nout-3]!=0) 
	    p_mc->Add(hlp=new VHAAG_ND(m_nin,m_nout,j,firsthaag));
	  if (!firsthaag) firsthaag=hlp;
	}
      }
      return 0;
    }

  };// end of class VHAAG_ND_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(VHAAG_ND_Channel_Generator,"VHAAG_ND",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,VHAAG_ND_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new VHAAG_ND_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    VHAAG_ND_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Vegas-improved HAAG integrator";
}
