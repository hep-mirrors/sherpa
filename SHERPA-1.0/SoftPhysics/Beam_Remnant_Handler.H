#ifndef Beam_Remnant_Handler_H
#define Beam_Remnant_Handler_H

#include "Remnant_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Blob_List.H"
#include "Primordial_KPerp.H"

namespace SHERPA {

  class Beam_Remnant_Handler {
  private:

    PDF::ISR_Handler           *p_isr;
    BEAM::Beam_Spectra_Handler *p_beam;

    SHERPA::Remnant_Base *p_beampart[2];
    Primordial_KPerp     *p_kperp;

    ATOOLS::Blob *p_beamblob[2];

    ATOOLS::Blob_List     *p_bloblist;
    ATOOLS::Particle_List *p_particlelist;

    std::string  m_path, m_file;
    bool         m_fill;

  public :

    // constructor
    Beam_Remnant_Handler(const std::string path,const std::string file,
			 PDF::ISR_Handler *const isr,
			 BEAM::Beam_Spectra_Handler *const beam,
			 const double scale);

    // destructor
    ~Beam_Remnant_Handler();
    
    // member functions
    bool FillBunchBlobs(ATOOLS::Blob_List *const bloblist,
			ATOOLS::Particle_List *const particlelist=NULL);
    bool FillBeamBlobs(ATOOLS::Blob_List *const bloblist,
		       ATOOLS::Particle_List *const particlelist=NULL);

    // inline functions
    inline Remnant_Base *const BeamParticle(const size_t i) const 
    { return p_beampart[i]; }

  };// end of class Beam_Remnant_Handler

}// end of namespace SHERPA

#endif

