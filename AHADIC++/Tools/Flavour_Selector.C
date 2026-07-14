#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include <cassert>

using namespace AHADIC;
using namespace ATOOLS;

Flavour_Selector::Flavour_Selector(Hadronisation_Reweighting * reweighting) :
  p_reweighting(reweighting) {}

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

  double disc {m_norms[0] * ran->Get()};
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

  // reweight with the selection probabilities of the different flavours,
  // including the different norms and popweights of the variations
  auto opt {m_options.find(ret)};
  if(opt == m_options.end())
    opt = m_options.find(ret.Bar());
  if(opt == m_options.end())
    THROW(fatal_error, "No flavour selected.");
  if(m_norms[0] == 0) return ret;
  p_reweighting->FlavourSelectionReweighting(opt->second->popweights, m_norms);
  p_reweighting->RecordFlavourPop(ret, Emax); // OUTPUT

  return ret;
}

void Flavour_Selector::Norm(const double & mmax,const bool & vetodi)
{
  std::fill(m_norms.begin(), m_norms.end(), 0);
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();fdit++) {
    if (vetodi && fdit->first.IsDiQuark()) continue;
    if (fdit->second->popweights[0]>0. && fdit->second->massmin<mmax/2.) {
      for (size_t ivar=0; ivar<m_n_variations; ++ivar)
	m_norms[ivar] += fdit->second->popweights[ivar];
    }
  }
}

void Flavour_Selector::Init() {
  Constituents * constituents(hadpars->GetConstituents());
  m_mmin = constituents->MinMass();
  m_mmax = constituents->MaxMass();
  m_mmin2 = ATOOLS::sqr(m_mmin);
  m_mmax2 = ATOOLS::sqr(m_mmax);
  m_n_variations = p_reweighting->NumberOfVariations();
  m_norms.resize(m_n_variations);
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
