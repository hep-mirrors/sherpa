#ifndef SHERPA_PerturbativePhysics_Hard_Decay_Handler_H
#define SHERPA_PerturbativePhysics_Hard_Decay_Handler_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Scoped_Settings.H"
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

  typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> ParticlePair;
  typedef std::pair<std::pair<ATOOLS::Particle *,ATOOLS::Particle *>,ATOOLS::Particle *> ParticlePairPair;
  typedef std::vector<ParticlePair> ParticlePair_Vector;
  typedef std::vector<ParticlePairPair> ParticlePairPair_Vector;

  class Hard_Decay_Handler : public Decay_Handler_Base {
    std::string              m_resultdir, m_offshell;
    std::set<std::string>    m_disabled_channels;
    std::map<ATOOLS::Flavour,std::set<std::string> > m_forced_channels;
    std::map<std::string, double> m_external_widths;
    int                      m_store_results;
    bool                     m_decay_tau, m_set_widths,
                             m_br_weights, m_usemass;
    double                   m_min_prop_width;
    ATOOLS::Flavour_Set      m_decmass;
    std::map<ATOOLS::Flavour, std::map<std::string, std::vector<double> > >
                             m_read;
    PDF::Cluster_Definitions_Base *p_clus;
    ATOOLS::NLO_subevtlist *p_newsublist;
    double                   m_int_accuracy;
    int                      m_int_niter, m_int_target_mode;

    void SetDecayMasses();
    void InitializeDirectDecays(PHASIC::Decay_Table* dt);
    void InitializeOffshellDecays(PHASIC::Decay_Table* dt);
    void SetHOSMWidths(ATOOLS::Scoped_Settings& s);
    bool TriggerOffshell(PHASIC::Decay_Channel* dc, std::vector<PHASIC::Decay_Channel*> new_dcs);
    std::vector<PHASIC::Decay_Channel*> ResolveDecay(PHASIC::Decay_Channel* dc);
    bool ProperVertex(MODEL::Single_Vertex* sv);
    
    void AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                            ATOOLS::Blob* blob,
                            size_t& imax,
                            size_t idmother);
    void AddSplitPhotonsClustering(ATOOLS::Cluster_Amplitude*& ampl,
                              const ATOOLS::Particle_Vector daughters,
                              ParticlePairPair_Vector& splitphotons,
                              size_t& imax,
                              const std::vector<size_t>& ids);
    void AddPhotonsClustering(ATOOLS::Cluster_Amplitude*& ampl,
                              const ATOOLS::Particle_Vector daughters,
                              ParticlePair_Vector& photons,
                              size_t& imax,
                              const std::vector<size_t>& ids);
    void UnsplitPhotons(const ATOOLS::Particle_Vector& splitphotonproducts,
                        ParticlePairPair_Vector& splitphotons);
    void AssignPhotons(const ATOOLS::Particle_Vector& daughters,
                       ParticlePair_Vector& photons);
    void AssignSplitPhotons(const ATOOLS::Particle_Vector& daughters,
                            ParticlePairPair_Vector& splitphotons);
    ATOOLS::Vec4D RecombinedMomentum(const ATOOLS::Particle * daughter,
                                     const ParticlePair_Vector& photons,
                                     const ParticlePairPair_Vector& splitphotons,
                                     size_t& stat);

    void ReadDecayTable(ATOOLS::Flavour decayer);
    void WriteDecayTables();
    bool CalculateWidth(PHASIC::Decay_Channel* dc);

    void FindDecayProducts(ATOOLS::Particle* decayer,
                           std::list<ATOOLS::Particle*>& decayproducts);
    double BRFactor(ATOOLS::Blob* blob) const;

  public :
    Hard_Decay_Handler();
    ~Hard_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Particle* inpart);
    void TreatInitialBlob(ATOOLS::Blob* blob,
                          METOOLS::Amplitude2_Tensor* amps,
                          const ATOOLS::Particle_Vector& origparts);

    bool DefineInitialConditions(ATOOLS::Cluster_Amplitude* ampl,
                                 ATOOLS::Blob* initial_blob);

    double Mass(const ATOOLS::Flavour &fl) const;
    bool Decays(const ATOOLS::Flavour& flav);

    inline void SetCluster(PDF::Cluster_Definitions_Base *clus) { p_clus=clus; }
  };
}

#endif
