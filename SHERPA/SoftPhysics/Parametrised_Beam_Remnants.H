#ifndef SHERPA_SoftPhysics_Parametrised_Beam_Remnants_H
#define SHERPA_SoftPhysics_Parametrised_Beam_Remnants_H

#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "SHERPA/SoftPhysics/Primordial_KPerp.H"
#include "ATOOLS/Org/Return_Value.H"

namespace SHERPA {

  class Parametrised_Beam_Remnants {
  private:
    PDF::ISR_Handler           * p_isr;
    PDF::Remnant_Base          * p_beampart[2];
    Primordial_KPerp           * p_kperp;

    ATOOLS::Blob               * p_beamblob[2];

    std::string                  m_path, m_file;
    short unsigned int           m_beam;
    void                         InitBeamBlob(const int beam);
  public :

    // constructor
    Parametrised_Beam_Remnants(const std::string path,const std::string file,
			       PDF::ISR_Handler *const isr,
			       BEAM::Beam_Spectra_Handler *const beam);

    // destructor
    ~Parametrised_Beam_Remnants();
    
    // member functions
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List *const,
					     ATOOLS::Particle_List *const=NULL);

    // inline functions
    inline PDF::Remnant_Base * BeamParticle(const size_t i) const { 
      return p_beampart[i]; 
    }
    inline void SetBeamBlob(ATOOLS::Blob * const blob,const int & beam) {
      p_beampart[beam]->Clear();
      p_beampart[beam]->ClearErrors();
      p_beamblob[beam] = blob;
    }
    inline void SetScale(const double & scale) {
      for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
    }
    inline void CleanUp() {
      for (short unsigned int i=0;i<2;++i) p_beampart[i]->Clear();
    }

  };// end of class Parametrised_Beam_Remnants

}// end of namespace SHERPA

#endif

