#ifndef SHERPA_SoftPhysics_Beam_Remnant_Handler_H
#define SHERPA_SoftPhysics_Beam_Remnant_Handler_H

#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "SHRiMPS/Main/Shrimps.H"

namespace SHERPA {
  class Beam_Remnant_Handler {
  private:
    REMNANTS::Remnant_Handler  * p_remnants, * p_bunchremnants;
    SHRIMPS::Shrimps           * p_shrimps;
    BEAM::Beam_Spectra_Handler * p_beam;
    ATOOLS::Blob               * p_beamblobs[2];
    Soft_Collision_Handler     * p_schandler;
    
    bool m_bunchrescatter, m_fill, m_vmode;
    int  m_beam;
    std::string m_name;

    ATOOLS::Return_Value::code TreatNoFill(ATOOLS::Blob_List *const bloblist);
    ATOOLS::Return_Value::code FillBunchBlobs(ATOOLS::Blob_List *const,
					      ATOOLS::Particle_List *const=NULL);
    bool           FillSimpleBunchBlobs(ATOOLS::Blob_List *const bloblist);
    bool           FillRescatterBunchBlobs(ATOOLS::Blob_List *const bloblist);
    ATOOLS::Blob * FillBunchBlob(int,ATOOLS::Particle *);
  public :
    Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler *const beam,
			 REMNANTS::Remnant_Handler *const remnants,
			 Soft_Collision_Handler *const softcollisions);
    ~Beam_Remnant_Handler();
    
    ATOOLS::Return_Value::code FillBeamAndBunchBlobs(ATOOLS::Blob_List *const,
						     const bool & onlyBunch=false);
    ATOOLS::Return_Value::code FillRescatterBeamBlobs(ATOOLS::Blob_List *const);
    void CleanUp(const size_t & mode=0);
    
    inline void AddBunchRescattering(REMNANTS::Remnant_Handler *const remnants,
				     Soft_Collision_Handler *const softcollisions) {
      if (remnants && remnants->Id()==PDF::isr::bunch_rescatter) {
	p_bunchremnants = remnants; p_schandler = softcollisions;
	m_bunchrescatter = true;
      }
    }
    inline void SetShrimps(SHRIMPS::Shrimps * shrimps) { 
      if (shrimps!=NULL) { p_shrimps = shrimps; m_name = std::string("Shrimps"); }
    }
    inline const bool NeedsToDealWithRescattering()  {
      return p_bunchremnants->NeedsToMakeBeamBlobs();
    }
    inline REMNANTS::Remnant_Handler * GetRemnants()      const { return p_remnants; }
    inline REMNANTS::Remnant_Handler * GetBunchRemnants() const { return p_bunchremnants; }
    inline const int Fill()                               const { return m_fill; }
    inline const std::string & Name()                     const { return m_name; }
   };
}

#endif

