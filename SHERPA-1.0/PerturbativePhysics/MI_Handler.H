#ifndef MI_Handler_H
#define MI_Handler_H

#include "Run_Parameter.H"
#include "Blob_List.H"

namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler;          }
namespace MODEL  { class Model_Base;           }
namespace AMISIC { class Amisic;               }

namespace SHERPA {

  class Matrix_Element_Handler;

  class MI_Handler {
  public:
    
    enum TypeID { None    =  0,
		  Amisic  =  1,
		  Unknown = 99 };

  private:

    AMISIC::Amisic   *p_amisic;

    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    TypeID      m_type;
    std::string m_name;

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
    
    SHERPA::Matrix_Element_Handler *HardMEHandler();
    SHERPA::Matrix_Element_Handler *SoftMEHandler();

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

    void SameHardProcess(ATOOLS::Blob *blob);
    void SameSoftProcess(ATOOLS::Blob *blob);

    bool GenerateHardEvent(ATOOLS::Blob_List *bloblist);
    bool GenerateSoftEvent(ATOOLS::Blob_List *bloblist);

    void SameHardEvent(ATOOLS::Blob_List *bloblist);
    void SameSoftEvent(ATOOLS::Blob_List *bloblist);

  }; // end of class MI_Handler

} // end of namespace SHERPA

#endif
