#include "Pair_Popper.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Pair_Popper::Pair_Popper(const double tension) :
  p_constituents(hadpars.GetConstituents()), 
  m_totweight(0),m_tension(tension),m_minmt2(1.e6)
{
  cout<<METHOD<<":"<<tension<<endl;
  double  mass2;
  Flavour flav;
  for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
       iter!=p_constituents->CCMap.end();iter++) { 
    m_totweight += iter->second->TotWeight();
    cout<<"   Add : "<<iter->first<<" "<<iter->second->TotWeight()<<endl;
    if (m_tension>0.) {
      flav  = iter->first;
      mass2 = sqr(p_constituents->Mass(flav));
      if (mass2<m_minmt2) m_minmt2 = mass2;
    }
  }
}


Return_Value::code Pair_Popper::Pop(Flavour & flav,double & mt2,const double upper)
{
  //cout<<METHOD<<":"<<m_totweight<<" for "<<upper<<","<<m_minmt2<<endl;
  double actweight(m_totweight),mass2,disc;
  flav = Flavour(kf::none);
  if (m_minmt2>upper/4.) return Return_Value::Error;
  do {
    do { mt2 = m_minmt2 - log(ran.Get())*m_tension/M_PI; } while (mt2>upper/4.);
    disc = actweight*ran.Get();
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) {
      flav = iter->first; 
      mass2 = sqr(p_constituents->Mass(flav));
      if (m_minmt2<=mass2 && mass2<mt2) disc -= iter->second->TotWeight();
      //cout<<" "<<mass2<<"("<<flav<<") < "<<m_minmt2<<" ? "<<(m_minmt2<=m_minmt2)
      //	  <<" & "<<mass2<<" <  mt2 = "<<mt2<<" -> disc : "<<disc<<endl;
      if (disc<=1.e-12) { 
	if (flav==Flavour(kf::c))
	  cout<<METHOD<<" "<<mass2<<"("<<flav<<") < "<<m_minmt2<<" ? "<<(m_minmt2<=m_minmt2)
	      <<" & "<<mass2<<" <  mt2 = "<<mt2<<" -> disc : "<<disc<<endl;
	if (flav.IsDiQuark()) flav = flav.Bar();
	return Return_Value::Success;
      }
    }
  } while (disc>0.);
  flav = (ran.Get()>0.5)?Flavour(kf::u):Flavour(kf::d);
  return Return_Value::Success;
}


Return_Value::code Pair_Popper::Pop(Flavour & flav,const double upper)
{
  double refmass(flav!=Flavour(kf::none)?min(upper,p_constituents->Mass(flav)):upper);
  double actweight = m_totweight;
  if (flav!=Flavour(kf::none)) {
    Flavour trial;
    actweight = 0.;
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) {
      trial = iter->first;
      if (p_constituents->Mass(trial)<refmass) actweight += iter->second->TotWeight();
    }
    if (actweight==0.) return Return_Value::Error;
  }
  double disc = actweight*ran.Get();
  for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
       iter!=p_constituents->CCMap.end();iter++) {
    disc -= iter->second->TotWeight();
    if (disc<=1.e-12) { 
      flav = iter->first; 
      if (flav.IsDiQuark()) flav = flav.Bar();
      return Return_Value::Success;
    }
  }
  msg.Error()<<"Warning in "<<METHOD<<" :"<<endl
	     <<"   Selecting did not work : tot = "<<m_totweight<<", disc = "<<disc<<"."<<endl
	     <<"   Will return u or d."<<endl;
  if (ran.Get()<=.5) flav = Flavour(kf::u);
                else flav = Flavour(kf::d);
  rvalue.IncWarning(METHOD);
  return Return_Value::Warning;
}
