#ifndef AHADIC_Tools_Flavour_Selector_H
#define AHADIC_Tools_Flavour_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  struct  DecaySpecs { double weight,massmin,popweight; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;

  class Flavour_Selector {
  private:
    FlavDecayMap m_options;
    double       m_mmin, m_mmin2, m_mmax, m_mmax2, m_sumwt;

    double Norm(const double & mmax=-1.,const bool & vetodi=false); 
  public:
    Flavour_Selector();
    ~Flavour_Selector();

    void   InitWeights();
    
    ATOOLS::Flavour operator()(const double & Emax,const bool & vetodi);
    const double MinimalMass() const { return m_mmin; }
  };
}

#endif
