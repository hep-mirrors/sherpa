#ifndef SHERPA_PerturbativePhysics_MI_Handler_H
#define SHERPA_PerturbativePhysics_MI_Handler_H

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace ATOOLS { class Cluster_Amplitude;    }
namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler;          }
namespace PHASIC { class Process_Base;         }
namespace MODEL  { class Model_Base;           }
#ifdef USING__Amisic
namespace AMISIC { class Amisic;               }
#endif

namespace SHERPA {

  class Matrix_Element_Handler;
  class Shower_Handler;

  class MI_Handler: public ATOOLS::Mass_Selector {
  public:
    
    enum TypeID { None    =  0,
		  Amisic  =  1,
		  Unknown = 99 };

  private:
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;
#ifdef USING__Amisic
    AMISIC::Amisic             *p_amisic;
#endif

    ATOOLS::Cluster_Amplitude *p_ampl;
    PHASIC::Process_Base      *p_proc;

    Shower_Handler *p_shower;

    TypeID      m_type;
    std::string m_name;

    double m_ycut; 

  public:

    // constructor
    MI_Handler(std::string path,std::string file,MODEL::Model_Base *model,
	       BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr);

    // destructor
    ~MI_Handler();
    
    // getters
    TypeID      Type();
    std::string Name();

    std::string MIGenerator();
    
    PDF::ISR_Handler *ISRHandler();

    // member functions
    void Reset();
    void CleanUp();

    void SetScaleMin(double scalemax,unsigned int i);
    void SetScaleMax(double scalemin,unsigned int i);

    double ScaleMin(unsigned int i);
    double ScaleMax(unsigned int i);

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    bool GenerateSoftProcess(ATOOLS::Blob *blob);
    bool GenerateEvent(ATOOLS::Blob_List *bloblist);

    bool VetoHardProcess(ATOOLS::Blob *blob);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration();

    double Mass(const ATOOLS::Flavour &fl) const;

    // inline functions
    inline double YCut() const { return m_ycut; }

    inline PHASIC::Process_Base *Process() const
    { return p_proc; }

    inline void SetShowerHandler(Shower_Handler *const sh) { p_shower=sh; }

    inline Shower_Handler *Shower() const { return p_shower; }

  }; // end of class MI_Handler

} // end of namespace SHERPA

#endif
