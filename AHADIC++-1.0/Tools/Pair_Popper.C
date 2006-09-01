#include "Pair_Popper.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Pair_Popper::Pair_Popper() :
  p_constituents(hadpars.GetConstituents()), 
  m_ptmode(2),m_flavmode(1),
  m_tension(hadpars.Get(string("Tension"))),
  m_ptexp(hadpars.Get(string("<pt/ptmax>"))),
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


Flavour Pair_Popper::SelectFlavour(const double upper,const bool nodiquarks)
{
  Flavour flav = Flavour(kf::none);
  if (m_minmass>upper) return flav;
  double mass, totweight(m_totweight);
  if (upper>m_maxmass) {
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) { 
      flav = iter->first;
      if (nodiquarks&&flav.IsDiQuark()) continue;
      mass = p_constituents->Mass(flav);
      if (2.*mass>upper) {
	switch (m_flavmode) {
	case 2:
	  // dynamic suppression
	  totweight -= iter->second->TotWeight() * exp(-M_PI*sqr(mass)/m_tension);
	  break;
	case 1:
	default:
	  // no dynamic suppression
	  totweight -= iter->second->TotWeight();
	  break;
	}
      }
    }
  }
  double wt(1.);
  while (wt>0.) {
    wt = totweight*ran.Get();
    for (FlavCCMap_Iterator iter=p_constituents->CCMap.begin();
	 iter!=p_constituents->CCMap.end();iter++) { 
      flav = iter->first;
      if (nodiquarks&&flav.IsDiQuark()) continue;
      mass = p_constituents->Mass(flav);
      if (2.*mass<upper) {
	switch (m_flavmode) {
	case 2:
	  wt -= iter->second->TotWeight()*exp(-M_PI*sqr(mass)/m_tension);
	  break;
	case 1:
	default:
	  wt -= iter->second->TotWeight();
	  break;
	}
      }
      if (wt<=1.e-12) { wt=-1.; break; }
    }
  }
  if (flav.IsDiQuark()) flav = flav.Bar();
  return flav;
}

double Pair_Popper::SelectPT(const double upper)
{
  double pt2;
  int    maxtrials(100);
  switch (m_ptmode) {
  case 2:
    // Shifted Gaussian
    while (maxtrials>0) {
      pt2 = -m_tension/M_PI*log(ran.Get());
      if (exp(-M_PI/m_tension*sqr(sqrt(pt2)-m_ptexp*sqrt(upper)))>ran.Get()) {
	if (pt2<upper) return sqrt(pt2);
      }
      maxtrials--;
    }
    break;
  case 1:
    // Gaussian
    while (maxtrials>0) {
      pt2 = -m_tension/M_PI*log(ran.Get());
      if (exp(-M_PI/m_tension*pt2)>ran.Get()) {
	if (pt2<upper) return sqrt(pt2);
      }
      maxtrials--;
    }
    break;
  case 0:
  default:
    // Exponential pt^2 ~ exp(-Pi/tension pt^2)
    while (maxtrials>0) {
      pt2 = -m_tension/M_PI*log(ran.Get());
      if (pt2<upper) return sqrt(pt2);
      maxtrials--;
    }
    break;
  }
  return sqrt(ran.Get()*upper);
}

double Pair_Popper::PopWeight(Flavour & flav) {
  return exp(-M_PI*sqr(p_constituents->Mass(flav))/m_tension);
}
