#ifndef Amisic_H
#define Amisic_H

#include "Simple_Chain.H"

namespace AMISIC {
  
  class Amisic: public File_IO_Base {
  public:

    enum ModelID { None        =  0,
		   // hard models
                   SimpleChain =  1,
		   // soft models (not yet implemented)
		   Unknown     = 99 };

  private:
    
    ModelID m_hardmodel, m_softmodel;
    MI_Base *p_hardbase, *p_softbase;

    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    bool m_external;

  public:

    // constructors
    Amisic();
    Amisic(MODEL::Model_Base *_p_model,
	   BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr);

    // destructor
    ~Amisic();

    // getters
    ModelID HardModel();
    ModelID SoftModel();

    EXTRAXS::XS_Base *HardXS();
    EXTRAXS::XS_Base *SoftXS();

    // member functions
    bool Initialize();
    void Reset();

    void SameHardProcess(ATOOLS::Blob *blob);
    void SameSoftProcess(ATOOLS::Blob *blob);

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    bool GenerateHardEvent(ATOOLS::Blob_List *blobs);
    
    bool GenerateSoftProcess(ATOOLS::Blob *blob);
    bool GenerateSoftEvent(ATOOLS::Blob_List *blobs);
    
    bool GenerateEvent(ATOOLS::Blob_List *blobs);

    bool SelectHardModel(ModelID _m_hardmodel);
    bool SelectSoftModel(ModelID _m_softmodel);

    static ModelID StringToModelID(std::string model);
    static std::string ModelIDToString(ModelID model);

  }; // end of class Amisic

} // end of namespace AMISIC

#include "Amisic.IC"

#endif
