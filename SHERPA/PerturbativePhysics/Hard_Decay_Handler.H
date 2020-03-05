#ifndef SHERPA_PerturbativePhysics_Hard_Decay_Handler_H
#define SHERPA_PerturbativePhysics_Hard_Decay_Handler_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHASIC++/Decays/Decay_Handler.H"
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
  class Decay_Table;
}

namespace SHERPA {

  class Hard_Decay_Handler : public PHASIC::Decay_Handler {
    std::string              m_resultdir, m_offshell;
    std::set<std::string>    m_disabled_channels;
    std::map<ATOOLS::Flavour,std::set<std::string> > m_forced_channels;
    std::map<std::string, double> m_external_widths;
    int                      m_store_results;
    bool                     m_set_widths,
                             m_br_weights, m_usemass;
    double                   m_min_prop_width;
    ATOOLS::Flavour_Set      m_decmass;
    std::map<ATOOLS::Flavour, std::map<std::string, std::vector<double> > >
                             m_read;
    ATOOLS::NLO_subevtlist *p_newsublist;
    double                   m_int_accuracy;
    int                      m_int_niter, m_int_target_mode;

    void SetDecayMasses();
    void InitializeDirectDecays(PHASIC::Decay_Table* dt);
    void InitializeOffshellDecays(PHASIC::Decay_Table* dt);
    void CustomizeDecayTables();
    bool TriggerOffshell(PHASIC::Decay_Channel* dc, std::vector<PHASIC::Decay_Channel*> new_dcs);
    std::vector<PHASIC::Decay_Channel*> ResolveDecay(PHASIC::Decay_Channel* dc);
    bool ProperVertex(MODEL::Single_Vertex* sv);
    
    void ReadDecayTable(ATOOLS::Flavour decayer);
    void WriteDecayTables();
    bool CalculateWidth(PHASIC::Decay_Channel* dc);

    void FindDecayProducts(ATOOLS::Particle* decayer,
                           std::list<ATOOLS::Particle*>& decayproducts);
    double BRFactor(ATOOLS::Blob* blob) const;

  public :
    Hard_Decay_Handler();
    ~Hard_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Blob_List* bloblist, ATOOLS::Particle* inpart);
    void AfterTreatInitialBlob(ATOOLS::Blob* blob);

    double Mass(const ATOOLS::Flavour &fl) const;
    inline std::string Name() { return "HardDecays"; }
  };
}

#endif
