#include "Basic_Func.H"
#include "Basic_Sfuncs.H"
#include "String_Generator.H"
#include "Model_Base.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Kabbala Basic_Pfunc::P(Pfunc* p1)
{ 
  if (p1->zerowidth==1) {
    p1->value = Pcalc(p1->fl,-1);
    return sgen->GetPnumber(p1,-1);
  }
  p1->value = Pcalc(p1->fl,p1->momnum);
  return sgen->GetPnumber(p1,p1->momnum);
}

Complex Basic_Pfunc::Pcalc(const Flavour& fl,const int a)
{ 
  if (a<0) return Complex(sqrt(M_PI/(fl.Mass()*fl.Width())),0.);
  return Propagator((BS->Momentum(a)).Abs2(),fl);
}

Complex Basic_Pfunc::Pcalc(const int fl,const int a)
{ 
  return Pcalc(Flavour((kf_code)(fl)),a);
}

Complex Basic_Pfunc::Propagator(double p2,Flavour fl)
{
  Complex value;
  
  if(fl.IsKK()){    
    if(MODEL::s_model->ScalarNumber(std::string("KK_mode"))>0) value=KKProp(p2);
    else {
      value = Complex(1.,0.)/
	Complex(p2-sqr(fl.Mass()),fl.Mass()*fl.Width());
    }
  }
  else {
    value = Complex(1.,0.)/
      Complex(p2-sqr(fl.Mass()),fl.Mass()*fl.Width());
  }
  //extra i
  if (fl.IsFermion() || fl.IsScalar() || fl.IsTensor()) value *= Complex (0.,1.);
  if (fl.IsVector())                                    value *= Complex (0.,-1.);

  return value;
}

double Basic_Pfunc::Ifunc(double x,int ed)
{
  double a=0.;
  if((ed%2)==0){
    for(int k=2;k<ed;k+=2)a-=pow(x,k)/k;
    return a-0.5*log(sqr(x)-1);
  }
  for(int k=1;k<ed;k+=2)a-=pow(x,k)/k;
  return a+0.5*log((x+1.)/(x-1.));
} 

double Basic_Pfunc::IEfunc(double x,int ed)
{
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

Complex Basic_Pfunc::KKProp(double p2)
{
  int    ed  = MODEL::s_model->ScalarNumber(std::string("ED"));
  double gn  = MODEL::s_model->ScalarConstant(std::string("G_Newton"));
  double ms  = MODEL::s_model->ScalarConstant(std::string("M_s"));
  double msq = MODEL::s_model->ScalarConstant(std::string("M2_s"));

  double vr,vv;
  switch(MODEL::s_model->ScalarNumber(std::string("KK_mode"))){
  case 1:
    if(ed==2)vr=log(msq/ATOOLS::dabs(p2));
    else vr=2./(ed-2);
    return Complex(-0.5*vr/sqr(msq)/gn,0.);

  case 2:
    if(p2>0){
      vr=ms/sqrt(p2);
      vv= 1./(pow(vr,ed+2)*sqr(p2)*gn);
      return Complex(Ifunc(vr,ed)*vv,-0.5*M_PI*vv);
    }
    vr=ms/sqrt(-p2);
    vv= 1./(pow(vr,ed+2)*sqr(p2)*gn);
    return Complex(-IEfunc(vr,ed)*vv,0.);    

  case 3:
    return Complex(-1./(M_PI*sqr(msq)*gn),0.);

  case 4:
    return Complex(1./(M_PI*sqr(msq)*gn),0.);
  case 5:
    return Complex(-.5/(sqr(msq)*gn),0.);
  }
  return Complex(0.,0.);
}



