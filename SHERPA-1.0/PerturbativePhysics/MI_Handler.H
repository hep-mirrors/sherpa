#ifndef MI_Handler_H
#define MI_Handler_H

#ifdef USING__Amisic
#define USING__SHERPA
#include "Amisic.H"
#else
#include "XS_Base.H"
#include "Model_Base.H"
#include "Run_Parameter.H"
#include "Blob_List.H"
#endif

namespace SHERPA {

  class Matrix_Element_Handler; 

}

namespace SHERPA {

  class MI_Handler {
  public:
    
    enum TypeID { None    =  0,
#ifdef USING__Amisic
		  Amisic  =  1,
#endif
		  Unknown = 99 };

  private:

#ifdef USING__Amisic
    AMISIC::Amisic   *p_amisic;
#endif

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

    PDF::ISR_Handler       *ISRHandler();

    // member functions
    void Reset();

    void SetScaleMin(double scalemax,unsigned int i);
    void SetScaleMax(double scalemin,unsigned int i);

    unsigned int NIn();
    unsigned int NOut();

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    bool GenerateSoftProcess(ATOOLS::Blob *blob);

    void SameHardProcess(ATOOLS::Blob *blob);
    void SameSoftProcess(ATOOLS::Blob *blob);

  }; // end of class MI_Handler

} // end of namespace SHERPA

#endif
