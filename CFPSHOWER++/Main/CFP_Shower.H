#ifndef CFPSHOWER_Main_CFP_Shower_H
#define CFPSHOWER_Main_CFP_Shower_H

#include "CFPSHOWER++/Shower/Shower.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace CFPSHOWER {
  class CFP_Shower : public PDF::Shower_Base {
  private:
    Shower                       m_shower;
    std::vector<Configuration *> m_configs;
    Cluster_Definitions        * p_cluster;
    ATOOLS::Mass_Selector      * p_massselector;

    double m_muQ2;
    bool   m_test;

    ATOOLS::Cluster_Amplitude * InitConfigs(ATOOLS::Cluster_Amplitude * const ampl);
    Configuration * Convert(ATOOLS::Cluster_Amplitude * const ampl,
			    std::map<ATOOLS::Cluster_Leg*,Parton*> & lmap);
    void ExtractParton(ATOOLS::Blob *const bl,Parton *const p);
    void TestShower(const double & E=45.6,
		    const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_u));
    Configuration * MakeSimpleFFConfiguration(const double & E=45.6,
					      const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_u));
  public:
    CFP_Shower(const PDF::Shower_Key &key);
    ~CFP_Shower();

    bool PrepareShower(ATOOLS::Cluster_Amplitude * const ampl,
		       const bool & soft=false);
    int  PerformShowers();
    int  PerformDecayShowers();
    bool ExtractPartons(ATOOLS::Blob_List *const bl);
    void CleanUp();

    PDF::Cluster_Definitions_Base * GetClusterDefinitions();
  };
}

#endif
