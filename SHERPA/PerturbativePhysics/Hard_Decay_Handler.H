#ifndef SHERPA_PerturbativePhysics_Hard_Decay_Handler_H
#define SHERPA_PerturbativePhysics_Hard_Decay_Handler_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"

namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *>  Vertex_List;
}

namespace ATOOLS {
  class NLO_subevtlist;
}

namespace PHASIC {
  class Decay_Channel;
}

namespace SHERPA {
  
  class Hard_Decay_Handler : public Decay_Handler_Base {
    std::string              m_path, m_file, m_resultdir, m_offshell;
    std::set<std::string>    m_disabled_channels;
    std::map<ATOOLS::Flavour,std::set<std::string> > m_forced_channels;
    bool                     m_store_results, m_decay_tau, m_set_widths, m_br_weights;
    std::map<ATOOLS::Flavour, std::map<std::string, std::vector<double> > > 
                             m_read;
    PDF::Cluster_Definitions_Base *p_clus;
    ATOOLS::NLO_subevtlist *p_newsublist;

    void InitializeDirectDecays(PHASIC::Decay_Table* dt);
    void RefineDecaysThreshold(PHASIC::Decay_Table* dt);
    void RefineDecaysByWidth(PHASIC::Decay_Table* dt);
    std::vector<PHASIC::Decay_Channel*> ResolveDecay(PHASIC::Decay_Channel* dc);
    bool ProperVertex(MODEL::Single_Vertex* sv);
    
    void AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                            ATOOLS::Blob* blob,
                            size_t& imax,
                            size_t idmother);
    
    void ReadDecayTable(ATOOLS::Flavour decayer);
    void WriteDecayTables();
    bool CalculateWidth(PHASIC::Decay_Channel* dc);

    void FindDecayProducts(ATOOLS::Particle* decayer,
                           std::list<ATOOLS::Particle*>& decayproducts);
    double BRFactor(ATOOLS::Blob* blob) const;

  public :
    Hard_Decay_Handler(std::string,std::string);
    ~Hard_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Particle* inpart);
    void TreatInitialBlob(ATOOLS::Blob* blob,
                          METOOLS::Amplitude2_Tensor* amps,
                          const ATOOLS::Particle_Vector& origparts);

    void DefineInitialConditions(ATOOLS::Cluster_Amplitude* ampl,
                                 ATOOLS::Blob* initial_blob);

    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.Mass(); }
    bool Decays(const ATOOLS::Flavour& flav);

    inline void SetCluster(PDF::Cluster_Definitions_Base *clus) { p_clus=clus; }
  };
}

#endif

