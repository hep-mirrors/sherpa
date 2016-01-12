#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Generator.H"
//#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;

using namespace std;

#define Complex_Mass_Scheme

Kabbala Basic_MassTermfunc::MassTerm(int a)
{
  if (iabs(a)<99) {
    cerr<<"Bad Prop in Mass_Term !!!"<<endl;
    abort();
  }

  Pfunc* p1 = NULL;
  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->arg[0]==iabs(a)) break;
  }

  assert(p1!=NULL);
  if ((p1->fl).Mass()==0.0) return Kabbala(string("1"),Complex(1.,0.));

  return sgen->GetMassnumber(Sign(a)*p1->momnum,p1->fl,MassTermCalc(Sign(a)*p1->momnum,p1->fl));
}

Complex Basic_MassTermfunc::MassTermCalc(int a,int fl)
{
  Flavour flav = Flavour((kf_code)(iabs(fl)));
  if (fl<0) flav = flav.Bar();
  return MassTermCalc(a,flav);
}

Complex Basic_MassTermfunc::MassTermCalc(int a,Flavour flav)
{
#ifdef Complex_Mass_Scheme
  Complex mass;
#else
  double mass;
#endif

  mass = flav.Mass();
  
#ifdef Complex_Mass_Scheme
  mass = sqrt(mass*mass-Complex(0.,1)*mass*flav.Width());
#endif

  if (a<0)           mass = -mass;
  if (flav.IsAnti()) mass = -mass;
  if (flav.MassSign()==-1) mass = -mass;

  return 1.+mass/csqrt((BS->Momentum(iabs(a))).Abs2());
}

