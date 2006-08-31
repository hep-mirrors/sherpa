#include "Pair_Popper.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Pair_Popper::Pair_Popper() :
  p_constituents(hadpars.GetConstituents()), 
  m_tension(hadpars.Get(string("Tension"))),
  m_mode(2),
  m_totweight(0), m_maxmass(0.), m_minmass(100.)
{
  double  wt(1.),mass;
  Flavour flav;
  for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
       iter!=p_constituents->CCMap.end();iter++) { 
    flav = iter->first;
    mass = p_constituents->Mass(flav);
    wt   = iter->second->TotWeight() * exp(-M_PI*sqr(mass)/m_tension);
    if (mass>m_maxmass) m_maxmass = mass;
    if (mass<m_minmass) m_minmass = mass;
    m_totweight += wt;
  }
}


Flavour Pair_Popper::SelectFlavour(const double upper)
{
  Flavour flav = Flavour(kf::none);
  if (m_minmass>upper) return flav;
  double mass, totweight(m_totweight);
  if (upper>m_maxmass) {
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) { 
      flav = iter->first;
      mass = p_constituents->Mass(flav);
      if (2.*mass>upper) {
	totweight -= iter->second->TotWeight() * exp(-M_PI*sqr(mass)/m_tension);
      }
    }
  }
  double wt(1.);
  while (wt>0.) {
    wt = totweight*ran.Get();
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) { 
      flav = iter->first;
      mass = p_constituents->Mass(flav);
      if (2.*mass<upper) {
	wt -= iter->second->TotWeight()*exp(-M_PI*sqr(mass)/m_tension);
      }
      if (wt<=1.e-12) { wt=-1.; break; }
    }
  }
  if (flav.IsDiQuark()) flav = flav.Bar();
  return flav;
}

double Pair_Popper::SelectPT(const double upper)
{
  double pt(1.1*upper), pt2;
  switch (m_mode) {
  case 2:
    // Shifted Gaussian
    while (pt>upper) {
      pt2 = sqrt(-m_tension/M_PI*log(ran.Get()));
      if (exp(-M_PI/m_tension*sqr(sqrt(pt2)-sqrt(0.25*upper)))>ran.Get()) pt = sqrt(pt2);
    }
    break;
  case 1:
    // Gaussian
    while (pt>upper) {
      pt2 = -m_tension/M_PI*log(ran.Get());
      if (exp(-M_PI/m_tension*pt2)>ran.Get()) pt = sqrt(pt2);
    }
    break;
  case 0:
  default:
    // Exponential pt^2 ~ exp(-Pi/tension pt^2)
    while (pt>upper) pt = sqrt(-m_tension/M_PI*log(ran.Get()));
    break;
  }
  return pt;
}

double Pair_Popper::PopWeight(Flavour & flav) {
  return exp(-M_PI*sqr(p_constituents->Mass(flav))/m_tension);
}
