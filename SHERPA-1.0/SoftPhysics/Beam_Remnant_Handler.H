#ifndef Beam_Remnant_Handler_H
#define Beam_Remnant_Handler_H

#include "Remnant_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Blob_List.H"
#include "Primordial_KPerp.H"
#include "Return_Value.H"

namespace SHERPA {

  class Matrix_Element_Handler;

  class Beam_Remnant_Handler {
  private:

    PDF::ISR_Handler           * p_isr;
    BEAM::Beam_Spectra_Handler * p_beam;

    PDF::Remnant_Base          * p_beampart[2];
    Primordial_KPerp           * p_kperp;

    ATOOLS::Blob               * p_beamblob[2];

    Matrix_Element_Handler     * p_mehandler;
    std::string                  m_path, m_file;
    bool                         m_fill;
    short unsigned int           m_beam;

    void                       Reset(ATOOLS::Blob_List *const);
    ATOOLS::Return_Value::code FillBunchBlobs(ATOOLS::Blob_List *const,
					      ATOOLS::Particle_List *const=NULL);
    ATOOLS::Blob *             FillBunchBlob(const int,ATOOLS::Particle *);
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List *const,
					     ATOOLS::Particle_List *const=NULL);
    bool                       FillBeamBlob(const int,ATOOLS::Particle *);

  public :

    // constructor
    Beam_Remnant_Handler(const std::string path,const std::string file,
			 PDF::ISR_Handler *const isr,
			 BEAM::Beam_Spectra_Handler *const beam);

    // destructor
    ~Beam_Remnant_Handler();
    
    // member functions
    ATOOLS::Return_Value::code FillBeamAndBunchBlobs(ATOOLS::Blob_List *const);
    void SetScale(const double scale);

    // inline functions
    inline PDF::Remnant_Base *const BeamParticle(const size_t i) const 
    { return p_beampart[i]; }

  };// end of class Beam_Remnant_Handler

}// end of namespace SHERPA

#endif

