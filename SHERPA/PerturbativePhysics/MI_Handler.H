#ifndef SHERPA_PerturbativePhysics_MI_Handler_H
#define SHERPA_PerturbativePhysics_MI_Handler_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace ATOOLS   { class Cluster_Amplitude;    }
namespace BEAM     { class Beam_Spectra_Handler; }
namespace PDF      { class ISR_Handler;          }
namespace REMNANTS { class Remnant_Handler;      }
namespace PHASIC   { class Process_Base;         }
namespace MODEL    { class Model_Base;           }
namespace AMISIC   { class Amisic;               }
namespace SHRIMPS  { class Shrimps;              }

namespace SHERPA {
  class Matrix_Element_Handler;
  class Shower_Handler;

  class MI_Handler {
  public:
    
    enum typeID { none    =  0,
		  amisic  =  1,
		  shrimps =  2,
		  unknown = 99 };

  private:
    PDF::ISR_Handler          * p_isr;
    REMNANTS::Remnant_Handler * p_remnants;
    AMISIC::Amisic            * p_amisic;
    SHRIMPS::Shrimps          * p_shrimps;
    ATOOLS::Cluster_Amplitude * p_ampl;
    PHASIC::Process_Base      * p_proc;
    Shower_Handler            * p_shower;

    bool        m_stop;
    typeID      m_type;
    std::string m_name;

    double m_ycut; 
    
    void InitAmisic(MODEL::Model_Base *model);
    void InitShrimps(MODEL::Model_Base *model);
  public:
    MI_Handler(MODEL::Model_Base *model, PDF::ISR_Handler *isr,
               REMNANTS::Remnant_Handler * remnant_handler);
    ~MI_Handler();
    
    bool InitialiseMPIs(const double & scale);
    void SetMaxEnergies(const double & E1,const double & E2);
    void ConnectColours(ATOOLS::Blob * showerblob);
    
    ATOOLS::Blob              * GenerateHardProcess();
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob * blob);
    void Reset();
    void CleanUp();
    
    const double ScaleMin() const;
    const double ScaleMax() const;
    const double ImpactParameter() const;
    const ATOOLS::Vec4D SelectPositionForScatter() const;

    inline const typeID      & Type() const { return m_type; }
    inline const std::string & Name() const { return m_name; }
    inline const bool        & Done() const { return m_stop; }
    inline const double      & YCut() const { return m_ycut; }
    const bool   IsMinBias() const;
    
    inline PDF::ISR_Handler          * ISRHandler() const { return p_isr; }
    inline REMNANTS::Remnant_Handler * Remnants()   const { return p_remnants; }
    inline PHASIC::Process_Base      * Process()    const { return p_proc; }
    inline Shower_Handler            * Shower()     const { return p_shower; }
    inline AMISIC::Amisic            * Amisic()     const { return p_amisic; }
    inline SHRIMPS::Shrimps          * Shrimps()    const { return p_shrimps; }
    
    inline void SetShowerHandler(Shower_Handler * const sh) { p_shower=sh; }
    inline void SetRemnantHandler(REMNANTS::Remnant_Handler *const rh) {
      p_remnants = rh;
    }


    bool   VetoScatter(ATOOLS::Blob *blob);
    void   SetMassMode(const int & massmode);
    int    ShiftMasses(ATOOLS::Cluster_Amplitude * ampl);
  }; 
}

#endif
