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
    run_mode::code    m_runmode;
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
    int  InitMinimumBiasEvent(ATOOLS::Blob_List * blobs);
    ATOOLS::Blob * GenerateEvent();

    void Test(const std::string & dirname);


    inline Colour_Generator * GetColourGenerator()  {
      return (p_inelastic!=NULL)?p_inelastic->GetColourGenerator():NULL;
    }
    inline Omega_ik * GetEikonal() const { return p_active->GetEikonal(); }
    inline const double & XSec()   const { return m_xsec; } 

    
    inline void SetMaxEnergies(const double & E1,const double & E2) {
      if (p_active==p_inelastic && p_inelastic!=NULL)
	p_inelastic->SetMaxEnergies(E1,E2);
    }
			   
    inline double B()               const { return p_active?p_active->B():-1.; }
    inline int    NLadders()        const { return p_active?p_active->NLadders():1; }
    inline bool   IsMinBias()       const { return p_active?p_active->IsMinBias():false; }
    inline bool   IsLastRescatter() const { return p_active?p_active->IsLastRescatter():false; }
    inline double Smin()            const { return p_active?p_active->Smin():-1.; }
    inline double MinKT2()          const { return p_active?p_active->KT2min():-1.; }
    inline double Yhat()            const { return p_active?p_active->Yhat():0.; }
    inline double Ymax()            const { return p_active?p_active->Ymax():0.; }
    inline double TMax()            const { return p_active?p_active->TMax():0.; }
  };
}
#endif
