#ifndef SINIC_Main_Event_Generator_H
#define SINIC_Main_Event_Generator_H

#include "SINIC++/Cross_Sections/Cross_Sections.H"
#include "SINIC++/Event_Generation/Elastic_Event_Generator.H"
#include "SINIC++/Event_Generation/Single_Diffractive_Event_Generator.H"
#include "SINIC++/Event_Generation/Double_Diffractive_Event_Generator.H"
#include "SINIC++/Event_Generation/Quasi_Elastic_Event_Generator.H"
#include "SINIC++/Event_Generation/Inelastic_Event_Generator.H"
#include "SINIC++/Beam_Remnants/Beam_Remnant_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"

namespace PDF    { class ISR_Handler; }

namespace SINIC {
  class Event_Generator {
  private:
    run_mode::code    m_runmode;
    weight_mode::code m_weightmode;

    Cross_Sections                     * p_cross;
    Elastic_Event_Generator            * p_elastic;
    Single_Diffractive_Event_Generator * p_sdiff;
    Double_Diffractive_Event_Generator * p_ddiff;
    Quasi_Elastic_Event_Generator      * p_qelastic;
    Inelastic_Event_Generator          * p_inelastic;
    Event_Generator_Base               * p_active;

    double m_xsec, m_minkt2;
  public:
    Event_Generator() : 
      p_elastic(NULL), p_sdiff(NULL), p_qelastic(NULL), p_inelastic(NULL) {};
    Event_Generator(const run_mode::code & runmode,
		    const weight_mode::code & weightmode);
    ~Event_Generator();

    void Initialise(Cross_Sections * cross,
		    Beam_Remnant_Handler * beams,
		    const int & test);
    bool DressShowerBlob(ATOOLS::Blob * blob);
    int  MinimumBiasEvent(ATOOLS::Blob_List * blobs);

    inline Omega_ik * GetEikonal() const { return p_active->GetEikonal(); }
    inline const double & XSec()   const { return m_xsec; } 
    inline const double   Smin()   const { return p_active->Smin(); }
    inline double ShowerMinKT2() {
      if (p_active) return ATOOLS::Max(p_active->Smin(),m_minkt2); 
      return m_minkt2;
    }
  };
}
#endif
