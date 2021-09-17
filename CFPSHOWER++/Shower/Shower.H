#ifndef CFPShower_Shower_Shower_H
#define CFPShower_Shower_Shower_H

#include "CFPSHOWER++/Shower/Kernels.H"
#include "CFPSHOWER++/Tools/Configuration.H"

namespace MODEL {
  class Model_Base;
  class Running_AlphaS;
  class Running_AlphaQED;
  class Single_Vertex;
}

namespace PDF {
  class PDF_Base;
  class ISR_Handler;
}

namespace CFPSHOWER {
  class Cluster_Definitions;

  class Shower {
  private:
    double   m_t0[2], m_t0min;
    int      m_kinscheme, m_kfactor, m_cplscheme, m_MEcorrs;
    long int m_nmax_em, m_nem;

    Cluster_Definitions         * p_cluster;
    const ATOOLS::Mass_Selector * p_msel;
    Configuration               * p_config;
    Splitting                   * p_winner;
    
    std::map<ATOOLS::Flavour, Kernels *> * p_kernels[5];
    std::vector<Splitting *>               m_splittings;
    std::vector<std::vector<double> >      m_integrals;
    std::vector<Kernels *>                 m_splitkernels;
    
    double m_weight, m_weightover;

    bool        Evolve(Parton * splitter);
    double      InitialiseIntegrals(Parton * splitter);
    Splitting * GenerateTestSplitting(Parton * splitter,const double & sum);
    bool        PerformSplitting();
    void        EstablishSpectators(Parton * parton);
    void        EstablishSoftPartners();
    void        AddWeight(const double & t);
  public:
    Shower();
    ~Shower();

    bool Init(MODEL::Model_Base * const,PDF::ISR_Handler * const);
    bool Evolve(Configuration * config);

    void Reset(); 
    inline void SetMassSelector(const ATOOLS::Mass_Selector * ms) { p_msel = ms; }
    inline const ATOOLS::Mass_Selector * GetMassSelector()     const { return p_msel; }
    inline Cluster_Definitions * const GetClusterDefinitions() const { return p_cluster; }
    inline const double & Weight() const { return m_weight; }
  };
}
#endif
