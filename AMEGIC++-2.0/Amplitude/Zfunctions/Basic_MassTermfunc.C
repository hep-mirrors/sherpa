#include "Basic_Func.H"
#include "String_Generator.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;

using namespace std;

#define Complex_Mass_Scheme

Kabbala Basic_MassTermfunc::MassTerm(int a)
{
  if (iabs(a)<99) {
    cerr<<"Bad Prop in Mass_Term !!!"<<endl;
    abort();
  }

  Pfunc* p1;
  for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->arg[0]==iabs(a)) break;
  }

  double mass = AORGTOOLS::rpa.consts.Mass(p1->fl,sqr(AORGTOOLS::rpa.gen.Ecms()));

  if (AMATOOLS::IsZero(mass)) 
#ifdef Kabbala_on
    return Kabbala(string("1"),Complex(1.,0.));
#else
    return Complex(1.,0.);
#endif    

  return sgen->Get_Massnumber(Sign(a)*p1->momnum,p1->fl,MassTermCalc(Sign(a)*p1->momnum,p1->fl));
}

Complex Basic_MassTermfunc::MassTermCalc(int a,int fl)
{
  Flavour flav = Flavour(kf::code(iabs(fl)));
  if (fl<0) flav = flav.bar();
  return MassTermCalc(a,flav);
}

Complex Basic_MassTermfunc::MassTermCalc(int a,Flavour flav)
{
#ifdef Complex_Mass_Scheme
  Complex mass;
#else
  double mass;
#endif

  mass = AORGTOOLS::rpa.consts.Mass(flav,sqr(AORGTOOLS::rpa.gen.Ecms()));

#ifdef Complex_Mass_Scheme
  mass -= Complex(0.,1./2.)*flav.width();
#endif

  if (a<0)           mass = -mass;
  if (flav.isanti()) mass = -mass;
  if (flav.get_mass_sign()==-1) mass = -mass;

  return 1.+mass/csqrt((BS->Momentum(iabs(a))).abs2());
}

