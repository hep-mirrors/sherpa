#include "Single_XS.H"
#include "Single_Process.H"
#include "Process_Group.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"
#include "Running_AlphaS.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>


using namespace EXTRAXS;
using namespace PHASIC;
using namespace AMEGIC;
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
  SetThres(sqr(_fl[2].Mass()+_fl[3].Mass()));
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

double Single_XS::operator()(AMATOOLS::Vec4D * _p) {
  double _s, _t, _u;
  _s = (_p[0]+_p[1]).Abs2();
  _t = (_p[0]-_p[2]).Abs2();
  _u = (_p[0]-_p[3]).Abs2();
  return this->operator()(_s, _t, _u); 
}

double Single_XS::operator()(double s,double t,double u) {
  AORGTOOLS::msg.Error()<<"Virtual Method : Single_XS::operator()."<<std::endl; 
  return 0.; 
}

void Single_XS::MakeBroker(ISR::ISR_Handler * isr, BEAM::Beam_Spectra_Handler * beam,
			   APHYTOOLS::Selector_Data * _seldata,AMEGIC::Process_Group * _broker) 
{
  Pol_Info * _plavs   = 0;
  Flavour *  _fl      = new Flavour[Nin()+Nout()];
  for (int i=0; i<(Nin()+Nout()); i++) _fl[i] = Flavs()[i];
  broker = new Single_Process(Nin(),Nout(),Flavs(),isr,beam,_seldata,2,
			      rpa.me.KFactorScheme(),rpa.me.ScaleScheme(),
			      _plavs, AMEGIC::XS_MODE);
  broker->SetXS(this);
  _broker->Add(broker);
  delete [] _fl;
}

