#ifndef AMISIC_Main_Amisic_H
#define AMISIC_Main_Amisic_H
/*!
  \file Amisic.H
  \brief Declares the class Amisic
*/

#include "AMISIC++/Main/MI_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

namespace AMISIC {
  /*!
    \namespace AMISIC
    \brief Houses all classes dedicated to the generation of 
    the underlying event.
    
    AMISIC++ is an accronym for <em>A</em> <em>M</em>ultiple 
    <em>I</em>nteraction <em>S</em>imulation <em>I</em>n
    <em>C++</em>. Accordingly the AMISIC namespace houses all 
    classes which are concerned with the generation of underlying 
    events.
  */
  
  class Amisic: public ATOOLS::File_IO_Base {
  private:
    
    std::string  m_hardmodel, m_softmodel;
    MI_Base     *p_hardbase, *p_softbase;

    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    bool m_external;

  public:

    // constructors
    Amisic();
    /*!
      \fn Amisic()
      \brief The default constructor. 
    */
    Amisic(MODEL::Model_Base *const model,
	   BEAM::Beam_Spectra_Handler *const beam,
	   PDF::ISR_Handler *const isr);
    /*!
      \fn Amisic(MODEL::Model_Base *const model,
      BEAM::Beam_Spectra_Handler *const beam,
      PDF::ISR_Handler *const isr)
      \brief The standard constructor.
      Creates the underlying event handler employing 
      external ISR, beam spectra and model handlers.
    */

    // destructor
    ~Amisic();

    // member functions
    bool Initialize();
    /*!
      \fn bool Initialize()
      \brief Initializes the underlying event handling.
    */
    void Reset();
    /*!
      \fn void Reset()
      \brief Resets the underlying event handler.
      This method is to be called before any event to be generated.
    */
    void CleanUp();
    /*!
      \fn void CleanUp()
      \brief Cleans the underlying event handler.
      This method is to be called after any event generated.
    */

    bool GenerateHardProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool GenerateHardProcess(ATOOLS::Blob *blob)
      \brief Generates one hard interaction, which is part of
      the hard underlying event.
    */
    bool GenerateHardEvent(ATOOLS::Blob_List *blobs);
    /*!
      \fn bool GenerateHardEvent(ATOOLS::Blob_List *blobs)
      \brief Generates a multitude of hard interactions
      which taken together form the hard underlying event.
    */
    void SameHardProcess(ATOOLS::Blob *blob);
    /*!
      \fn void SameHardProcess(ATOOLS::Blob *blob)
      \brief Returns last hard intreraction.
    */
    
    bool GenerateSoftProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool GenerateSoftProcess(ATOOLS::Blob *blob)
      \brief Generates one soft interaction, which is part of
      the soft underlying event.
    */
    bool GenerateSoftEvent(ATOOLS::Blob_List *blobs);
    /*!
      \fn bool GenerateSoftEvent(ATOOLS::Blob_List *blobs)
      \brief Generates a multitude of soft interactions
      which taken together form the soft underlying event.
    */
    void SameSoftProcess(ATOOLS::Blob *blob);
    /*!
      \fn void SameSoftProcess(ATOOLS::Blob *blob)
      \brief Returns last soft intreraction.
    */
    
    bool VetoHardProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool VetoHardProcess(ATOOLS::Blob *blob)
      \brief Vetoes the hard interaction according to the 
      probability density employed in the hard underlying event model.
    */

    bool GenerateEvent(ATOOLS::Blob_List *blobs);
    /*!
      \fn bool GenerateEvent(ATOOLS::Blob_List *blobs)
      \brief Generates one underlying event, using the active
      hard and soft underlying event models.
    */

    bool SelectHardModel(const std::string &hardmodel);
    /*!
      \fn bool SelectHardModel(const std::string &hardmodel)
      \brief Selects the hard underlying event model 
      according to the given tag.
    */
    bool SelectSoftModel(const std::string &softmodel);
    /*!
      \fn bool SelectSoftModel(const std::string &softmodel)
      \brief Selects the soft underlying event model
      according to the given tag.
    */

    // inline functions
    inline const std::string &HardModel() const { return m_hardmodel; }
    inline const std::string &SoftModel() const { return m_softmodel; }

    inline MI_Base *HardBase() const { return p_hardbase; }
    inline MI_Base *SoftBase() const { return p_softbase; }

    inline PHASIC::Process_Base *HardXS() { return p_hardbase->XS(); }
    inline PHASIC::Process_Base *SoftXS() { return p_softbase->XS(); }

  }; // end of class Amisic
  /*!
    \class Amisic
    \brief The overall steering of the AMISIC package
    
    This class implements the overall steering of the 
    underlying event generation. It is interfaced by the Sherpa
    framework but it may also act as a stand-alone underlying event 
    generator. It updates the MI_Base instances with the information 
    which is necessary for the next event generation and delivers 
    the created underlying event as a blob to the outside world.
  */

} // end of namespace AMISIC

#endif
