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

Single_XS::Single_XS(int _nin,int _nout,Flavour * _fl)
{
  Init(_nin,_nout,_fl);
  GenerateName();
}

Single_XS::~Single_XS() {
  if (broker) delete broker;
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
    name    += string(fl[i].Name());
    if ((fl[i].Kfcode()==kf::e)   ||
	(fl[i].Kfcode()==kf::mu)  ||
	(fl[i].Kfcode()==kf::tau) ||
	(fl[i].Kfcode()==kf::Hmin)) {
      //kill last
      name.erase(name.length()-1,1);
      if (fl[i].IsAnti()) name += string("+");
                     else name += string("-");      
    }
    else {
      if (fl[i].IsAnti()) name += string("b"); 
    }
    name += string("_");
  }
  name.erase(name.length()-1,1);

  name      += string(" -> ");
  for (int i=nin;i<nin+nout;i++) {
    name    += string(fl[i].Name());
    if ((fl[i].Kfcode()==kf::e)   ||
	(fl[i].Kfcode()==kf::mu)  ||
	(fl[i].Kfcode()==kf::tau) ||
	(fl[i].Kfcode()==kf::Hmin)) {
      //kill last
      name.erase(name.length()-1,1);
      if (fl[i].IsAnti()) name += string("+");
                       else name += string("-");      
    }
    else {
      if (fl[i].IsAnti()) name += string("b"); 
    }
    name += string("_");
  }
  name.erase(name.length()-1,1);
}

double Single_XS::operator()(double s,double t,double u) {
  AORGTOOLS::msg.Error()<<"Virtual Method : Single_XS::operator()."<<std::endl; 
  return 0.; 
}
