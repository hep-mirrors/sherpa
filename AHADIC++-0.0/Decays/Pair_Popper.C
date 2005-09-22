#include "Pair_Popper.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Pair_Popper::Pair_Popper() :
  p_constituents(hadpars.GetConstituents()), m_totweight(0)
{
  for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
       iter!=p_constituents->CCMap.end();iter++) { 
    //cout<<"Pop "<<iter->first<<" "<<iter->second->TotWeight()<<endl;
    m_totweight += iter->second->TotWeight();
  }
  //abort();
}


bool Pair_Popper::Pop(Flavour & flav)
{
  double refmass = 100.,   actweight = m_totweight;
  if (flav!=Flavour(kf::none)) {
    Flavour trial;
    actweight = 0.;
    refmass   = p_constituents->Mass(flav);
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) {
      trial = iter->first;
      if (p_constituents->Mass(trial)<refmass) actweight += iter->second->TotWeight();
    }
    if (actweight==0.) return false;
  }
  double disc = actweight*ran.Get();
  for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
       iter!=p_constituents->CCMap.end();iter++) {
    disc -= iter->second->TotWeight();
    if (disc<=1.e-12) { 
      flav = iter->first; 
      if (flav.IsDiQuark()) flav = flav.Bar();
      return true; 
    }
  }
  msg.Error()<<"Warning in Pair_Popper::Pop :"<<endl
	     <<"   Selecting did not work : tot = "<<m_totweight<<", disc = "<<disc<<"."<<endl
	     <<"   Will return u or d."<<endl;
  if (ran.Get()<=.5) flav = Flavour(kf::u);
                else flav = Flavour(kf::d);
  return true; 
}
