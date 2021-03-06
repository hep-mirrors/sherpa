#ifndef SHRIMPS_Main_Event_Generator_H
#define SHRIMPS_Main_Event_Generator_H

#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
//#include "SHRiMPS/Event_Generation/Elastic_Event_Generator.H"
//#include "SHRiMPS/Event_Generation/Single_Diffractive_Event_Generator.H"
//#include "SHRiMPS/Event_Generation/Double_Diffractive_Event_Generator.H"
//#include "SHRiMPS/Event_Generation/Quasi_Elastic_Event_Generator.H"
#include "SHRiMPS/Event_Generation/Inelastic_Event_Generator.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"

namespace PDF    { class ISR_Handler; }

namespace SHRIMPS {
  class Cluster_Algorithm;
  
  class Event_Generator {
  private:
    run_mode::code    m_runmode, m_thisevent;
    weight_mode::code m_weightmode;

    Inelastic_Event_Generator * p_inelastic;
    Event_Generator_Base      * p_active;

    double m_xsec, m_minkt2;
  public:
    Event_Generator();
    ~Event_Generator();

    void Reset();
    void Initialise();
    bool DressShowerBlob(ATOOLS::Blob * blob);
    int  MinimumBiasEvent(ATOOLS::Blob_List * blobs);

    void Test(const std::string & dirname);


    inline Omega_ik * GetEikonal() const { return p_active->GetEikonal(); }
    inline const double & XSec()   const { return m_xsec; } 
    inline double         Smin()   const { return p_active->Smin(); }

    inline void SetCluster(Cluster_Algorithm * cluster) {
      if (p_inelastic) p_inelastic->SetCluster(cluster);
    }

    inline double ShowerMinKT2() {
      if (p_active) return p_active->Smin()*m_minkt2; 
      return 0.;
    }
    inline bool IsLastRescatter() const {
      if (p_active) return p_active->IsLastRescatter();
      return false;
    }
    inline double TMax() const {
      if (p_active) return p_active->TMax();
      return 0.;
    }
    inline int NLadders() const {
      if (p_active) return p_active->NLadders();
      return 1;
    }
  };
}
#endif
