#include "Single_XS.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"
#include "Running_AlphaS.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>


using namespace EXTRAXS;
using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

int fak(int N)
{
  if (N == 0) return 1;
  if (N < 0) return 0;  
  int res =1;
  for (int i=1;i<=N;i++) res *= i;
  return res;
}

/*-------------------------------------------------------------------------------

  Constructor, generic sequence of flavours, name generation,
  symmetry factors (flavour, helicities ..)

  ------------------------------------------------------------------------------- */

Single_XS::Single_XS(int _nin,int _nout,Flavour * _fl)
{
  Init(_nin,_nout,_fl);
  GenerateName();
}

void Single_XS::GenerateName() {
  char help[20];
  sprintf(help,"%i",nin);
  name       = string(help);
  name      += string("_");
  sprintf(help,"%i",nout);
  name      += string(help);
  name      += string("_");

  for (int i=0;i<nin;i++) {
    name    += string(fl[i].name());
    if ((fl[i].kfcode()==kf::e)   ||
	(fl[i].kfcode()==kf::mu)  ||
	(fl[i].kfcode()==kf::tau) ||
	(fl[i].kfcode()==kf::Hmin)) {
      //kill last
      name.erase(name.length()-1,1);
      if (fl[i].isanti()) name += string("+");
                     else name += string("-");      
    }
    else {
      if (fl[i].isanti()) name += string("b"); 
    }
    name += string("_");
  }
  name.erase(name.length()-1,1);

  name      += string(" -> ");
  for (int i=nin;i<nin+nout;i++) {
    name    += string(fl[i].name());
    if ((fl[i].kfcode()==kf::e)   ||
	(fl[i].kfcode()==kf::mu)  ||
	(fl[i].kfcode()==kf::tau) ||
	(fl[i].kfcode()==kf::Hmin)) {
      //kill last
      name.erase(name.length()-1,1);
      if (fl[i].isanti()) name += string("+");
                       else name += string("-");      
    }
    else {
      if (fl[i].isanti()) name += string("b"); 
    }
    name += string("_");
  }
  name.erase(name.length()-1,1);
}

double Single_XS::Differential(vec4d * p) {
  s = (p[0]+p[1]).abs2();
  t = (p[0]-p[2]).abs2();
  u = (p[0]-p[3]).abs2();
  return DSigma(s,t,u);
};

double Single_XS::Differential2() {
  return DSigma2();
};

double Single_XS::DSigma(double s,double t,double u)
{
  lastdxs  = operator()(s,t,u);
  lastlumi = isr->Weight(fl);
  return last = lastdxs * lastlumi;
}

double Single_XS::DSigma2() { 
  if (fl[0]==fl[1]) return 0.;
  double tmp = lastdxs * isr->Weight2(fl); 
  last      += tmp;
  return tmp;
}

double Single_XS::operator()(double s,double t,double u) {
  AORGTOOLS::msg.Error()<<"Virtual Method : Single_XS::operator()."<<std::endl; 
  return 0.; 
}

void Single_XS::SetISR(ISR::ISR_Handler * _isr) { 
  msg.Debugging()<<"Single_XS::SetISR("<<_isr->Type()<<") for "<<name<<" : "<<_isr<<endl;
  isr = _isr; 
}

void Single_XS::SetBeam(BEAM::Beam_Handler * _beam) { 
  msg.Debugging()<<"Single_XS::SetISR("<<_beam->Type()<<") for "<<name<<" : "<<_beam<<endl;
  beam = _beam; 
}
