#include "Basic_Func.H"
#include "String_Generator.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Kabbala Basic_Pfunc::P(Pfunc* p1)
{ 
  p1->value = Pcalc(p1->fl,p1->momnum);
  return sgen->Get_Pnumber(p1,p1->momnum);
}

Kabbala Basic_Pfunc::P(const int& pntemp)
{
#ifdef Kabbala_on
  Kabbala value(string("1"),Complex(1.,0.));
#else
  Complex value(1.,0.);
#endif    
  //testcase
  //return value;
  
  /*cout<<" ++++++++++++++++++++++++++++++++++++++++ "<<endl;
  for (short int i=0;i<pntemp;i++) {
    cout<<i<<" : "<<ps[i].numb<<"   ";
    }*/

  for (short int i=0;i<pntemp;i++) {
    Pfunc* p1;
    int hit = 0;
    for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
      p1 = *pit;
      if (p1->momnum==ps[i].numb && (p1->fl).kfcode()==ps[i].kfcode) {
	hit = 1;
	break;
      }
    }
    if (hit) {
      if (p1->arg[0]>99 && !p1->fl.isscalar() && p1->on==0) {
	//cout<<"Multipliy inner with "<<p1->arg[0]<<" = "<<p1->fl<<endl;
	value *= P(p1);
      }
    }
  }

  return value;
}

Complex Basic_Pfunc::Pcalc(const Flavour& fl,const int& a)
{ return Propagator((BS->Momentum(a)).abs2(),fl);}

Complex Basic_Pfunc::Pcalc(const int& fl,const int& a)
{ return Pcalc(Flavour(kf::code(fl)),a);}

Complex Basic_Pfunc::Propagator(double p2,Flavour fl)
{
  Complex value = Complex(1.,0.)/
    Complex(p2-sqr(AORGTOOLS::rpa.consts.Mass(fl,sqr(AORGTOOLS::rpa.gen.Ecms()))),
	    AORGTOOLS::rpa.consts.Mass(fl,sqr(AORGTOOLS::rpa.gen.Ecms()))*
	    AORGTOOLS::rpa.consts.Width(fl,sqr(AORGTOOLS::rpa.gen.Ecms())));

  //extra i
  if (fl.isfermion() || fl.isscalar()) value *= Complex (0.,1.);
  if (fl.isvector())                   value *= Complex (0.,-1.);

  return value;
}







