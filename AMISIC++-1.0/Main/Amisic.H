#ifndef Amisic_H
#define Amisic_H

#include "Simple_Chain.H"
#include "XS_Base.H"

namespace AMISIC {
  
  class Amisic: public ATOOLS::File_IO_Base {
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
    Amisic(MODEL::Model_Base *const model,
	   BEAM::Beam_Spectra_Handler *const beam,
	   PDF::ISR_Handler *const isr);

    // destructor
    ~Amisic();

    // member functions
    bool Initialize();
    void Reset();
    void CleanUp();

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    bool GenerateHardEvent(ATOOLS::Blob_List *blobs);
    void SameHardProcess(ATOOLS::Blob *blob);
    void SameHardEvent(ATOOLS::Blob_List *blobs);
    
    bool GenerateSoftProcess(ATOOLS::Blob *blob);
    bool GenerateSoftEvent(ATOOLS::Blob_List *blobs);
    void SameSoftProcess(ATOOLS::Blob *blob);
    void SameSoftEvent(ATOOLS::Blob_List *blobs);
    
    bool GenerateEvent(ATOOLS::Blob_List *blobs);

    bool SelectHardModel(ModelID _m_hardmodel);
    bool SelectSoftModel(ModelID _m_softmodel);

    static ModelID StringToModelID(std::string model);
    static std::string ModelIDToString(ModelID model);

    // inline functions
    inline ModelID HardModel() const { return m_hardmodel; }
    inline ModelID SoftModel() const { return m_softmodel; }

    inline MI_Base *const HardBase() const { return p_hardbase; }
    inline MI_Base *const SoftBase() const { return p_softbase; }

#ifdef USING__Sherpa
    inline SHERPA::Matrix_Element_Handler *const HardMEHandler()
    { return p_hardbase->MEHandler(); }
    inline SHERPA::Matrix_Element_Handler *const SoftMEHandler()
    { return p_softbase->MEHandler(); }
#endif

    inline EXTRAXS::XS_Base *const HardXS() { return p_hardbase->XS(); }
    inline EXTRAXS::XS_Base *const SoftXS() { return p_softbase->XS(); }

  }; // end of class Amisic

} // end of namespace AMISIC

#endif
