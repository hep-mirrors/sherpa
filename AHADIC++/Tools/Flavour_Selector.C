#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
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
  ATOOLS::Flavour ret;

  // update norms
  Norm(Emax,vetodi);

  double disc {norms[0] * ran->Get()};
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (vetodi && fdit->first.IsDiQuark()) continue;
    if (fdit->second->popweights[0]>0. && fdit->second->massmin<Emax/2.)
      disc -= fdit->second->popweights[0];
    if (disc<=0.) {
      // have to bar flavours for diquarks
      ret = fdit->first.IsDiQuark()?fdit->first.Bar():fdit->first;
      break;
    }
  }

  // compute probabilities for different flavours here
  // will include different norms and popweights
  auto opt {m_options.find(ret)};
  if(opt == m_options.end())
    opt = m_options.find(ret.Bar());
  if(opt == m_options.end())
    THROW(fatal_error, "No flavour selected.");
  if(norms[0] == 0) return ret;
  const double p0 {opt->second->popweights[0] / norms[0]};
  std::cout << "DEBUG: FLAVS: " << ret << " "<< int(ret) <<" (";
  for(int i; i<opt->second->popweights.size(); ++i) {
    tmp_variation_weights[i] *= (opt->second->popweights[i] / norms[i]) / p0;
    std::cout << (opt->second->popweights[i] / norms[i]) / p0;
    if(i < opt->second->popweights.size() -1)
      std::cout << ",";
  }
  std::cout << ")" << std::endl;

  if(int(ret) == -3303) {
    std::cout << "DEBUG_ss" << std::endl;
    for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
      std::cout << fdit->first << std::endl;
    }
    std::cout << norms << std::endl;
    std::cout << ret << std::endl;
    std::cout << opt->first << std::endl;
    for(int i; i<opt->second->popweights.size(); ++i) {

      std::cout << "wgt = " << opt->second->popweights[i] / norms[i] << std::endl;
      std::cout << "p0  = " << p0 << std::endl;
      std::cout << "frac = " << tmp_variation_weights[i] << std::endl;
    }
  }

  return ret;
}

void Flavour_Selector::Norm(const double & mmax,const bool & vetodi)
{
  std::fill(norms.begin(), norms.end(), 0);
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (vetodi && fdit->first.IsDiQuark()) continue;
    if (fdit->second->popweights[0]>0. && fdit->second->massmin<mmax/2.) {
      for(int i{0}; i<norms.size(); ++i)
	norms[i] += fdit->second->popweights[i];
    }
  }
}

void Flavour_Selector::Init() {
  Constituents * constituents(hadpars->GetConstituents());
  m_mmin = constituents->MinMass();
  m_mmax = constituents->MaxMass();
  m_mmin2 = ATOOLS::sqr(m_mmin);
  m_mmax2 = ATOOLS::sqr(m_mmax);
  norms.resize(constituents->m_nvars);
  variation_weights.resize(constituents->m_nvars);
  DecaySpecs * decspec;
  for (FlavCCMap_Iterator fdit=constituents->CCMap.begin();
       fdit!=constituents->CCMap.end();fdit++) {
    if (!fdit->first.IsAnti()) {
      decspec = new DecaySpecs;
      decspec->popweight  = constituents->TotWeight(fdit->first);
      decspec->massmin    = constituents->Mass(fdit->first);
      decspec->popweights = constituents->Weights(fdit->first);
      m_options[fdit->first] = decspec;
    }
  }
}
