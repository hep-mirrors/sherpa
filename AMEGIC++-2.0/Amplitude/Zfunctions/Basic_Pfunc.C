#include "Basic_Func.H"
#include "String_Generator.H"
#include "Run_Parameter.H"
#include "Couplings_LED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Kabbala Basic_Pfunc::P(Pfunc* p1)
{ 
  p1->value = Pcalc(p1->fl,p1->momnum);
  return sgen->Get_Pnumber(p1,p1->momnum);
}


Complex Basic_Pfunc::Pcalc(const Flavour& fl,const int& a)
{ return Propagator((BS->Momentum(a)).Abs2(),fl);}

Complex Basic_Pfunc::Pcalc(const int& fl,const int& a)
{ return Pcalc(Flavour(kf::code(fl)),a);}

Complex Basic_Pfunc::Propagator(double p2,Flavour fl)
{
  Complex value;
  if(fl.IsKK()){
    //Model dependent, to be improved!!!
    Couplings_LED  CplLED;
    CplLED.Init();
    if(CplLED.DoSum()) value=CplLED.KKProp(p2);
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







