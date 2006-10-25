#ifndef MI_Handler_H
#define MI_Handler_H

#include "Run_Parameter.H"
#include "Blob_List.H"
#include "CXXFLAGS.H"

namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler;          }
namespace MODEL  { class Model_Base;           }
#ifdef USING__Amisic
namespace AMISIC { class Amisic;               }
#endif

namespace SHERPA {

  class Matrix_Element_Handler;

  class MI_Handler {
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
    Matrix_Element_Handler     *p_hardmehandler, *p_softmehandler;

    TypeID      m_type;
    std::string m_name;

    int    m_scalescheme;
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
    
    Matrix_Element_Handler *HardMEHandler() const 
    { return p_hardmehandler; }
    Matrix_Element_Handler *SoftMEHandler() const 
    { return p_softmehandler; }

    PDF::ISR_Handler *ISRHandler();

    // member functions
    void Reset();
    void CleanUp();

    void SetScaleMin(double scalemax,unsigned int i);
    void SetScaleMax(double scalemin,unsigned int i);

    unsigned int NIn();
    unsigned int NOut();

    double ScaleMin(unsigned int i);
    double ScaleMax(unsigned int i);

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    bool GenerateSoftProcess(ATOOLS::Blob *blob);
    bool GenerateEvent(ATOOLS::Blob_List *bloblist);

    bool VetoHardProcess(ATOOLS::Blob *blob);
    
    // inline functions
    inline int ScaleScheme() const { return m_scalescheme; }

    inline double YCut() const { return m_ycut; }

  }; // end of class MI_Handler

} // end of namespace SHERPA

#endif
