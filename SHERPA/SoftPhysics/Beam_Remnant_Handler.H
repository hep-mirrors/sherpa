#ifndef SHERPA_SoftPhysics_Beam_Remnant_Handler_H
#define SHERPA_SoftPhysics_Beam_Remnant_Handler_H

#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
//#include "SHRiMPS/Main/Shrimps.H"

namespace SHERPA {
  class Beam_Remnant_Handler {
  private:
    REMNANTS::Remnant_Handler  * p_remnants;
    //SHRIMPS::Shrimps           * p_shrimps;
    BEAM::Beam_Spectra_Handler * p_beam;
    ATOOLS::Blob               * p_beamblobs[2];

    int m_beam, m_fill, m_vmode;
    std::string m_name;

    ATOOLS::Return_Value::code TreatNoFill(ATOOLS::Blob_List *const bloblist);
    ATOOLS::Blob             * FillBunchBlob(const int,ATOOLS::Particle *);
    ATOOLS::Return_Value::code FillBunchBlobs(ATOOLS::Blob_List *const,
					      ATOOLS::Particle_List *const=NULL);
  public :
    Beam_Remnant_Handler(const std::string path,const std::string file,
			 BEAM::Beam_Spectra_Handler *const beam,
			 REMNANTS::Remnant_Handler *const remnants,
			 Soft_Collision_Handler *const softcollisions);
    ~Beam_Remnant_Handler();
    ATOOLS::Return_Value::code FillBeamAndBunchBlobs(ATOOLS::Blob_List *const);
    void CleanUp(const size_t & mode=0);

    inline const int Fill() const { return m_fill; }
    inline const std::string & Name() const { return m_name; }
  };
}

#endif

