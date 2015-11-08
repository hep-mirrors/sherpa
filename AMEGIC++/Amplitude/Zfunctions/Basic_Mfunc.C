#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/String/String_Generator.H"
//#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Kabbala Basic_Mfunc::M(const int a)
{
  Pfunc* p1 = 0;
  
  int hit = 0;

  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->momnum==ps[iabs(a)].numb && (p1->fl).Kfcode()==ps[iabs(a)].kfcode) {
      hit = 1;
      break;
    }
  }
  if (hit==0) return sgen->GetEnumber(0.);

  Complex mass2 = Complex(0.,0.);

  //double mass = 0.;
  
  if (p1->arg[0]>99) {
      mass2 = Complex(sqr(p1->fl.Mass()),0);
      if (p1->fl.Width()>0.) 
	  mass2 -= Complex(0,p1->fl.Mass()*
	  p1->fl.Width());
  }
  if (mass2==0.0) return sgen->GetEnumber(0.);

  return sgen->GetMnumber(p1->fl,1./mass2);

}
