#include "Basic_Func.H"
#include "String_Generator.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace AMATOOLS;
using namespace std;

Kabbala Basic_Mfunc::M(const int &a)
{
  Pfunc* p1 = 0;
  
  int hit = 0;

  for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->momnum==ps[iabs(a)].numb && (p1->fl).kfcode()==ps[iabs(a)].kfcode) {
      hit = 1;
      break;
    }
  }
  if (hit==0) return sgen->Get_Enumber(0.);

  Complex mass2 = Complex(0.,0.);

  //double mass = 0.;
  
  if (p1->arg[0]>99) {
      mass2 = Complex(sqr(AORGTOOLS::rpa.consts.Mass(p1->fl,sqr(AORGTOOLS::rpa.gen.Ecms()))),0);
      if (p1->fl.width()>0.) 
	  mass2 -= Complex(0,AORGTOOLS::rpa.consts.Mass(p1->fl,sqr(AORGTOOLS::rpa.gen.Ecms()))*
	  p1->fl.width());
  }
  if (AMATOOLS::IsZero(mass2)) return sgen->Get_Enumber(0.);

  return sgen->Get_Enumber(1./mass2);

}
