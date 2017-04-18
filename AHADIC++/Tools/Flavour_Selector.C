#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include <cassert>

using namespace AHADIC;
using namespace ATOOLS;

Flavour_Selector::Flavour_Selector() {}

Flavour_Selector::~Flavour_Selector() {
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) 
    delete fdit->second;
  m_options.clear();
}
  
ATOOLS::Flavour Flavour_Selector::
operator()(const double & Emax,const bool & vetodi) {
  double disc = Norm(Emax,vetodi) * ran->Get();
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (fdit->second->popweight>0. && fdit->second->massmin<Emax/2. &&
	!(vetodi && fdit->first.IsDiQuark())) 
      disc -= fdit->second->popweight;
    if (disc<=0.) {
      // have to bar flavours for diquarks
      return fdit->first.IsDiQuark()?fdit->first.Bar():fdit->first;
    }
  }
}

double Flavour_Selector::Norm(const double & mmax,const bool & vetodi) 
{
  double sumwt(0.), wt;
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (vetodi && fdit->first.IsDiQuark()) continue;
    if (fdit->second->popweight>0. && fdit->second->massmin<mmax/2.) {
      wt = fdit->second->popweight; 
      sumwt += wt;
    }   
  } 
  return sumwt;
}

void Flavour_Selector::InitWeights() {
  Constituents * constituents(hadpars->GetConstituents());
  m_mmin = constituents->MinMass();
  m_mmax = constituents->MaxMass();
  m_mmin2 = ATOOLS::sqr(m_mmin);
  m_mmax2 = ATOOLS::sqr(m_mmax);
  DecaySpecs * decspec;
  for (FlavCCMap_Iterator fdit=constituents->CCMap.begin();
       fdit!=constituents->CCMap.end();fdit++) {
    if (!fdit->first.IsAnti()) {
      decspec = new DecaySpecs;
      decspec->popweight = constituents->TotWeight(fdit->first);
      decspec->massmin   = constituents->Mass(fdit->first);
      m_options[fdit->first] = decspec;
    }
  }
  m_sumwt = Norm();
}
