#ifndef SHERPA_SoftPhysics_Beam_Remnant_Handler_H
#define SHERPA_SoftPhysics_Beam_Remnant_Handler_H

#include "PDF/Remnant/Remnant_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "SHERPA/SoftPhysics/Primordial_KPerp.H"
#include "ATOOLS/Org/Return_Value.H"

namespace SHERPA {

  class Beam_Remnant_Handler {
  private:

    PDF::ISR_Handler           * p_isr;
    BEAM::Beam_Spectra_Handler * p_beam;

    PDF::Remnant_Base          * p_beampart[2];
    Primordial_KPerp           * p_kperp;

    ATOOLS::Blob               * p_beamblob[2];

    std::string                  m_path, m_file;
    bool                         m_fill;
    short unsigned int           m_beam;
    int                          m_vmode, m_on;

    ATOOLS::Return_Value::code FillBunchBlobs(ATOOLS::Blob_List *const,
					      ATOOLS::Particle_List *const=NULL);
    ATOOLS::Blob *             FillBunchBlob(const int,ATOOLS::Particle *);
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List *const,
					     ATOOLS::Particle_List *const=NULL);
    void                       InitBeamBlob(const int beam);

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
    
    void CleanUp();
    
    // inline functions
    inline PDF::Remnant_Base *const BeamParticle(const size_t i) const 
    { return p_beampart[i]; }

  };// end of class Beam_Remnant_Handler

}// end of namespace SHERPA

#endif

