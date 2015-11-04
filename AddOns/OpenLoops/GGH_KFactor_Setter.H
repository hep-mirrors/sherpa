#ifndef PHASIC_GGH_KFactor_H
#define PHASIC_GGH_KFactor_H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "AddOns/OpenLoops/GGH_Process_Manager.H"

namespace ATOOLS {
  class NLO_subevt;
}

namespace AMEGIC{
  class Single_Process;
  class Single_Process_External;
}

namespace PHASIC {

  class GGH_KFactor_Setter : public KFactor_Setter_Base {
    
    typedef std::vector<Subprocess_Info> SubprocInfo_Vector;
    
  public:

    GGH_KFactor_Setter(const KFactor_Setter_Arguments &args);
    virtual ~GGH_KFactor_Setter();
    
    // Default KFactor method
    double KFactor(const int mode=0);
    // KFactor for Comix, not yet tested or validated
    double KFactor(const ATOOLS::NLO_subevt& evt);

  private:

    // Vertex correction for gg->H with
    // no kinematic dependency
    bool               m_set_vcorrection;
    double             vertex_correction;

    // is p_proc* a real correction?
    bool                m_real_corr;

    // Minimum pt threshold for massless partons:
    // Cluster back to lower multipl process below
    static const double s_ir_co;

    // default processes for the ratio tree_ME/loop_ME
    AMEGIC::Single_Process*          p_default_tree;
    AMEGIC::Single_Process_External* p_default_loop;

    double ClusterMassCorrectionFactor(ATOOLS::NLO_subevtlist* rsprocs);
    double ClusterMassCorrectionFactor();
    bool   IsCollinear(const ATOOLS::Vec4D_Vector& mom) const;
    void   PrintAuxProcs();
    void   SetOSVertexCorrection();
    double OSVertexCorrection();

  protected:

    // Default amplitude representing partonic process
    ATOOLS::Cluster_Amplitude* p_ampl;
    // 'Clustered' amplitude used for collinear configs
    ATOOLS::Cluster_Amplitude* p_next_ampl;

    // default reweighting factor
    virtual double        MassCorrectionFactor(const ATOOLS::Vec4D_Vector& p);
    // reweighting factor for collinear configs
    virtual double        MassCorrectionFactor(const ATOOLS::Cluster_Amplitude& ampl);

    // Reweighting factors for Comix, not yet tested and validated
    double                MassCorrectionFactor(const ATOOLS::NLO_subevt& evt);
    double                MassCorrectionFactor(const std::string& name,
					       const ATOOLS::Vec4D_Vector& p);
    // Used for setting default amplitude p_ampl
    virtual               ATOOLS::Cluster_Amplitude* GetAmpl() const;
    // Set p_next_ampl to a clustered configuration
    virtual void          SetNextAmplitude();
    
    virtual               ATOOLS::Vec4D_Vector GetMomenta() const;
    bool                  ContainsDecays(const PHASIC::Process_Base& proc) const;

  };

  extern GGH_Process_Manager s_procmanager;
}

namespace OpenLoops{
  class OpenLoops_Interface;
  extern OpenLoops_Interface* s_interface;
}

#endif
