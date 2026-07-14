#ifndef AHADIC_Tools_Flavour_Selector_H
#define AHADIC_Tools_Flavour_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  class Hadronisation_Reweighting;

  struct  DecaySpecs {
    double weight,massmin,popweight;
    std::vector<double> popweights;
  };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;

  class Flavour_Selector {
  private:
    FlavDecayMap m_options;
    double       m_mmin, m_mmin2, m_mmax, m_mmax2;
    Hadronisation_Reweighting * p_reweighting;
    size_t m_n_variations {0};
    std::vector<double> m_norms;

    void Norm(const double & mmax=-1.,const bool & vetodi=false);
  public:
    Flavour_Selector(Hadronisation_Reweighting * reweighting);
    ~Flavour_Selector();

    void   Init();

    ATOOLS::Flavour operator()(const double & Emax,const bool & vetodi);
    const double MinimalMass() const { return m_mmin; }
  };
}

#endif
