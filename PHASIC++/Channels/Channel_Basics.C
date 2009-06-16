#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

void Channel_Basics::Rotat(int lflag,Vec4D& p1 ,Vec4D p2,double** rot)
{
  if (lflag==0) {
    short int i,k,l;
    double r[2][3][3],pm[2],sp[2],cp[2],st[2],ct[2];
    Vec4D pp[2];
    pm[0] = Vec3D(p1).Abs();
    pm[1] = Vec3D(p2).Abs();
    pp[0] = (1./pm[0])*p1;
    pp[1] = (1./pm[1])*p2;
    for (i=0;i<2;i++) {
      ct[i] = pp[i][3];
      st[i] = sqrt(1.-sqr(ct[i]));
      if (IsEqual(dabs(ct[i]),1.)) {
	cp[i] = 1.;
	sp[i] = 0.;
      }
      else {
	cp[i] = pp[i][2]/st[i];
	sp[i] = pp[i][1]/st[i];
      }
      r[i][0][0]=  cp[i]; 
      r[i][0][1]=  sp[i]*ct[i]; 
      r[i][0][2]=  st[i]*sp[i]; 
      r[i][1][0]=  -sp[i]; 
      r[i][1][1]=  ct[i]*cp[i];   
      r[i][1][2]=  cp[i]*st[i];
      r[i][2][0]=  0.;
      r[i][2][1]=  -st[i]; 
      r[i][2][2]=  ct[i];
    }
    for (i=0;i<3;i++) {
      for (l=0;l<3;l++) {
	rot[i][l] = 0.;
	for (k=0;k<3;k++) 
	  rot[i][l] += r[0][i][k]*r[1][l][k];
      }
    }
    
    Vec4D p1new;
    p1new[0] = p2[0];
    for (short int i=0;i<3;i++) {
      p1new[i+1] = 0.;
      for (short int j=0;j<3;j++) {
	p1new[i+1] += rot[i][j]*p2[j+1];
      }
    }
  }
  else {
    short int i,j;
    p1[0] = p2[0];
    for (i=0;i<3;i++) {
      p1[i+1] = 0.;
      for (j=0;j<3;j++) {
	p1[i+1] += rot[i][j]*p2[j+1];
      }
    }
  }
}

void Channel_Basics::Boost(int lflag,Vec4D q,Vec4D& ph,Vec4D& p)
{
  if (q.Abs2() < 0.) {
    msg_Error()<<"Channel_Basics::Boost : Spacelike four vector ..."<<endl;
    return;
  }
  double rsq = sqrt(q.Abs2());
  if (lflag==0) {
    p[0] = (q[0]*ph[0]+Vec3D(q)*Vec3D(ph))/rsq;
    double c1 = (ph[0]+p[0])/(rsq+q[0]);
    p = Vec4D(p[0],Vec3D(ph)+c1*Vec3D(q));  
  }
  else {
    ph[0] = q*p/rsq;
    double c1 = (p[0]+ph[0])/(rsq+q[0]);
    ph = Vec4D(ph[0],Vec3D(p)-c1*Vec3D(q));  
  }
}

double Channel_Basics::SqLam(double s,double s1, double s2)
{
  double arg1 = sqr(s-s1-s2)-4.*s1*s2;
  if (arg1>0.) return sqrt(arg1)/s;
  msg_Error()<<"Channel_Basics::SqLam argument "<<arg1<<" <0 in Channel_Basics::sqlam()"<<endl
			<<"s;s1;s2: "<<s<<";"<<s1<<";"<<s2<<endl;
  return 0.;
}

double Channel_Basics::PeakedDist(double a,double cn,double cxm,double cxp,int k,double ran)
{
  double ce  = 1.-cn;
  double res = 0.;
  if (!IsZero(ce)) {
    res = k * (pow(ran*pow(a+k*cxp,ce)+(1.-ran)*pow(a+k*cxm,ce),1/ce)-a);
  }
  else {
//     if (cxp<-a) res = -k * (exp(ran*log(-k*a-cxp)+(1.-ran)*log(-k*a-cxm))-a);
//            else res = k * (exp(ran*log(a+k*cxp)+(1.-ran)*log(a+k*cxm))-a);
    res = k *( (a+k*cxm)*pow( (a+k*cxp)/(a+k*cxm) , ran) - a);
  }
  return res;
}

double Channel_Basics::PeakedGrid(double a,double cn,double cxm,double cxp,double res,int k,double &ran)
{
  double ce  = 1.-cn;
  if (!IsZero(ce)) {
    double amin = pow(a+k*cxm,ce);
    ran = (pow(a+k*res,ce)-amin)/(pow(a+k*cxp,ce)-amin);
  }
  else {
    double amin = a+k*cxm;
    ran = log((a+k*res)/amin)/log((a+k*cxp)/amin);
  }
  return res;
}

double Channel_Basics::PeakedWeight(double a,double cn,
				    double cxm,double cxp,int k)
{
  double ce = 1.-cn;
  double wt;
  if (!IsZero(ce)) wt = (pow(a+k*cxp,ce)-pow(a+k*cxm,ce))/(k*ce);
              else wt = log((a+k*cxp)/(a+k*cxm))/k;
  return wt;
}

double Channel_Basics::PeakedWeight(double a,double cn,
				    double cxm,double cxp,double res,int k,double &ran)
{
  double ce = 1.-cn;
  double wt;
  if (!IsZero(ce)) {
    double amin = pow(a+k*cxm,ce);
    wt = pow(a+k*cxp,ce)-amin;
    ran = (pow(a+k*res,ce)-amin)/wt;
    wt /= k*ce;
  }
  else {
    double amin = a+k*cxm;
    wt = log((a+k*cxp)/amin);
    ran = log((a+k*res)/amin)/wt;
    wt /= k;
  }
  return wt;
}

double Channel_Basics::BoundaryPeakedDist(double cxm,double cxp,double ran)
  //  1/(x(1-x))
{ 
  double fxp=1./cxp-1.;
  double fxm=1./cxm-1.;
  double pw = pow(fxm/fxp,ran);
  return pw/(fxm+pw);
}

double Channel_Basics::BoundaryPeakedWeight(double cxm,double cxp,double res,double &ran)
  //  1/(x(1-x))
{
  double fxp=1./cxp-1.;
  double fxm=1./cxm-1.;
  double wt=log(fxm/fxp);
  ran = log(fxm/(1./res-1.))/wt;
  return wt;
}

double Channel_Basics::Tj1(double cn,double amcxm,double amcxp,double ran)
{
  double ce= 1.-cn;
  double res = 0.;
  if (!ATOOLS::IsZero(ce)) res = pow(ran*pow(amcxm,ce)+(1.-ran)*pow(amcxp,ce),1./ce);
  else {
    if(amcxp>0.) res =  exp(ran*log(amcxm)+(1.-ran)*log(amcxp));
            else res = -exp(ran*log(-amcxm)+(1.-ran)*log(-amcxp));
  }
  return res;
}

double Channel_Basics::Hj1(double cn,double amcxm,double amcxp)
{
  double ce= 1.-cn;
  if (!ATOOLS::IsZero(ce)) return (pow(amcxp,ce)-pow(amcxm,ce))/ce;
  return log(amcxp/amcxm);
}



double Channel_Basics::PseudoAngleCut(double m1_sq,double E1,
				      double m2_sq,double E2) 
{
  double mu1_sq = m1_sq/sqr(E1);
  double mu2_sq = m2_sq/sqr(E2);
  double beta1  = sqrt(1.-mu1_sq);
  double beta2  = sqrt(1.-mu2_sq);
  double del1   = (mu1_sq-mu2_sq)/(beta1+beta2);
  double del2   = 1.-E2/E1;
  double arg;
  if (del2<0.1) {
    arg = 0.;
    for (short int i=1;i<56;i++) arg += pow(del2,i+1);
  }
  else arg= 1./(1.-del2)-1.-del2;
  return (del1*(del1*E2/E1 - 2.*beta1*del2)-mu1_sq*arg)/2.;
  //checked
}

double Channel_Basics::
FlatDist(double alpha,double min,double max,double R)
{
  double p=1.0+alpha;
  double Imin=pow(log(min),p), Imax=pow(log(max),p);
  if (min<1.0) return exp(-pow(Imax*R+(1.0-R)*Imin,1.0/p));
  return exp(pow(Imax*R+(1.0-R)*Imin,1.0/p));
}

double Channel_Basics::
FlatGrid(double alpha,double min,double max,double s,double &R)
{
  double p=1.0+alpha;
  double Imin=pow(log(min),p), Imax=pow(log(max),p);
  R=(pow(log(s),p)-Imin)/(Imax-Imin);
  return s;
}

double Channel_Basics::FlatWeight(double alpha,double min,double max)
{
  double p=1.0+alpha;
  double Imin=pow(log(min),p), Imax=pow(log(max),p);
  return (Imax-Imin)/p;
}

double Channel_Basics::
FlatWeight(double alpha,double min,double max,double s,double &R)
{
  double p=1.0+alpha;
  double Imin=pow(log(min),p), Imax=pow(log(max),p);
  R=(pow(log(s),p)-Imin)/(Imax-Imin);
  return (Imax-Imin)/p;
}


