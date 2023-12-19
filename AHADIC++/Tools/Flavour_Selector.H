#ifndef AHADIC_Tools_Flavour_Selector_H
#define AHADIC_Tools_Flavour_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
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
    std::vector<double> variation_weights;
    std::vector<double> tmp_variation_weights;
    std::vector<double> norms;

    void Norm(const double & mmax=-1.,const bool & vetodi=false);
  public:
    Flavour_Selector();
    ~Flavour_Selector();

    void   Init();

    void   accept_splitting() {
      for(int i{0}; i<variation_weights.size(); ++i)
	variation_weights[i] *= tmp_variation_weights[i];
    }

    std::vector<double>& get_variationweights() {
      return variation_weights;
    }

    void reset_var_weights() {
      std::fill(tmp_variation_weights.begin(), tmp_variation_weights.end(), 1);
    }

    void reset_variationweights(int nvars) {
      tmp_variation_weights.resize(nvars);
      variation_weights.resize(nvars);
      std::fill(variation_weights.begin(), variation_weights.end(), 1);
    };

    ATOOLS::Flavour operator()(const double & Emax,const bool & vetodi);
    const double MinimalMass() const { return m_mmin; }
  };
}

#endif
