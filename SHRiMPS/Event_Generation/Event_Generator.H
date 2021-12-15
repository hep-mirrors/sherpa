#ifndef SHRIMPS_Main_Event_Generator_H
#define SHRIMPS_Main_Event_Generator_H

#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Event_Generation/Elastic_Event_Generator.H"
#include "SHRiMPS/Event_Generation/Soft_Diffractive_Event_Generator.H"
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

    Inelastic_Event_Generator        * p_inelastic;
    Elastic_Event_Generator          * p_elastic;
    Soft_Diffractive_Event_Generator * p_soft_diffractive;
    Event_Generator_Base             * p_active;

    double m_xsec, m_xsec_inel, m_xsec_elas, m_xsec_diff;
    bool   m_eventsuccessful, m_mustinit;

    void InitGenerator(Cross_Sections * xsecs,const bool & test);
  public:
    Event_Generator(Cross_Sections * xsecs,const bool & test);
    ~Event_Generator();

    void Initialise(Remnant_Handler * remnants,Cluster_Algorithm * cluster);
    void Reset();
    bool DressShowerBlob(ATOOLS::Blob * blob);
    int  MinimumBiasEvent(ATOOLS::Blob_List * blobs);

    void Test(const std::string & dirname);


    inline Colour_Generator * GetColourGenerator()  {
      return (p_inelastic!=NULL)?p_inelastic->GetColourGenerator():NULL;
    }
    inline Omega_ik * GetEikonal() const { return p_active->GetEikonal(); }
    inline const double & XSec()   const { return m_xsec; } 
    inline double         Smin()   const { return p_active->Smin(); }
    inline double         B()      const { return p_active->B(); }
    
			   
    inline double MinKT2() {
      if (p_active) return p_active->KT2min();
      return 0.;
    }
    inline double Ymax() {
      if (p_active) return p_active->Ymax();
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
