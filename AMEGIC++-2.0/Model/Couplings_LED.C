#include <math.h>
#include <fstream> 
#include "Couplings_LED.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;



void Couplings_LED::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  ed            = dr.GetValue<int>("ED");
  gn            = dr.GetValue<double>("G_Newton");
  ms            = dr.GetValue<double>("M_s");
  rad           = dr.GetValue<double>("Radius");
  kk_sum_mode   = dr.GetValue<int>("KK-Sum");

  if(IsZero(ms)&&rad>0.){
    double gam;
    if(ed%2==0)gam=1.;
    else gam=sqrt(M_PI);
    for(int i=2-ed%2;i<ed;i+=2)gam*=0.5*i;

    ms=pow(gam*pow(4.*M_PI,.5*ed)/pow(rad,1.*ed)/gn,1./(2.+ed));
  }
}

double Couplings_LED::Gn(){return gn;}
double Couplings_LED::Ms(){return ms;}
double Couplings_LED::Msq(){return sqr(ms);}
int Couplings_LED::Ned(){return ed;}
double Couplings_LED::Kappa(){return sqrt(8.*M_PI*gn);}
double Couplings_LED::Omega(){return sqrt(4.*(-1.+ed)/(3.*(2.+ed)));}
double Couplings_LED::R()
{
  double gam;
  if(ed%2==0)gam=1.;
  else gam=sqrt(M_PI);
  for(int i=2-ed%2;i<ed;i+=2)gam*=0.5*i;

  return pow(gam*pow(4.*M_PI,.5*ed)/pow(ms,2.+(double(ed)))/gn,1./(double(ed)));
}

double Couplings_LED::Ifunc(double x)
{
  double a=0.;
  if((ed%2)==0){
    for(int k=2;k<ed;k+=2)a-=pow(x,k)/k;
    return a-0.5*log(sqr(x)-1);
  }
  for(int k=1;k<ed;k+=2)a-=pow(x,k)/k;
  return a+0.5*log((x+1.)/(x-1.));
} 

double Couplings_LED::IEfunc(double x)
{
  int n;
  double a=0.;
  if((ed%2)==0){
    for(int k=2;k<ed;k+=2){
      if(((k/2)%2)==0) a+=pow(x,k)/k;
      else             a-=pow(x,k)/k;
    }
    a+=0.5*log(sqr(x)+1);
    if((ed%4)==2) return a;
    else          return -a;
  }

  for(int k=1;k<ed;k+=2){
    if((((k+1)/2)%2)==0) a+=pow(x,k)/k;
    else                 a-=pow(x,k)/k;
  }
  a+=atan(x);
  if((ed%4)==1) return a;
  else          return -a;
} 

Complex Couplings_LED::KKProp(double p2)
{
  double vr,vv;
  switch(kk_sum_mode){
  case 1:
    if(ed==2)vr=log(Msq()/AMATOOLS::dabs(p2));
    else vr=2./(ed-2);
    return Complex(-0.5*vr/sqr(Msq())/gn,0.);

  case 2:
    if(p2>0){
      vr=ms/sqrt(p2);
      vv= 1./(pow(vr,ed+2)*sqr(p2)*gn);
      return Complex(Ifunc(vr)*vv,-0.5*M_PI*vv);
    }
    vr=ms/sqrt(-p2);
    vv= 1./(pow(vr,ed+2)*sqr(p2)*gn);
    return Complex(-Ifunc(vr)*vv,0.);    

  case 3:
    return Complex(-1./(M_PI*sqr(Msq())*gn),0.);

  case 4:
    return Complex(1./(M_PI*sqr(Msq())*gn),0.);
  case 5:
    return Complex(-.5/(sqr(Msq())*gn),0.);
  }
}


