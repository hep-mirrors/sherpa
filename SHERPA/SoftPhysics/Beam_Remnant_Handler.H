#ifndef SHERPA_SoftPhysics_Beam_Remnant_Handler_H
#define SHERPA_SoftPhysics_Beam_Remnant_Handler_H

#include "SHERPA/SoftPhysics/Parametrised_Beam_Remnants.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "SHRiMPS/Main/Shrimps.H"

namespace SHERPA {
  class Beam_Remnant_Handler {
  private:
    Parametrised_Beam_Remnants * p_parametrised;
    SHRIMPS::Shrimps               * p_shrimps;
    BEAM::Beam_Spectra_Handler * p_beam;
    ATOOLS::Blob               * p_beamblobs[2];

    int m_beam, m_fill, m_vmode;

    ATOOLS::Blob             * InitBeamBlob(const int beam);
    ATOOLS::Blob             * FillBunchBlob(const int,ATOOLS::Particle *);
    ATOOLS::Return_Value::code FillBunchBlobs(ATOOLS::Blob_List *const,
					      ATOOLS::Particle_List *const=NULL);
  public :
    Beam_Remnant_Handler(const std::string path,const std::string file,
			 BEAM::Beam_Spectra_Handler *const beam,
			 PDF::ISR_Handler *const isr,
			 Soft_Collision_Handler *const softcollisions);
    ~Beam_Remnant_Handler();
    ATOOLS::Return_Value::code FillBeamAndBunchBlobs(ATOOLS::Blob_List *const);
    void CleanUp(const size_t & mode=0);

    inline const int Fill() const { return m_fill; }
  };// end of class Beam_Remnant_Handler

}// end of namespace SHERPA

#endif

