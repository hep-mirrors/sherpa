#ifndef REMNANTS_Main_Parametrised_H
#define REMNANTS_Main_Parametrised_H

#include "REMNANTS/Main/Remnant_Handler.H"
#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace BEAM {
  class Beam_Spectra_Handler;
}

namespace REMNANTS {
  class Remnant_Handler;
}

namespace REMNANTS {
  class Parametrised {
  private:
    Remnant_Handler  * p_rhandler;
    ATOOLS::Blob     * p_beamblobs[2];
  public :
    Parametrised(const std::string path,const std::string file,
		 REMNANTS::Remnant_Handler * const remnants,
		 BEAM::Beam_Spectra_Handler *const beam);
    ~Parametrised();
    
    inline Remnant_Base * BeamParticle(const size_t i) const { 
      return p_rhandler->GetRemnant(i); 
    }
    inline void SetScale2(const double & scale2=-1.0) {
      p_rhandler->SetScale2(scale2);
    }
    inline void SetBeamBlob(ATOOLS::Blob * blob,const unsigned int & beam) {
      if (beam<2) {
	p_beamblobs[beam] = blob;
	p_rhandler->SetBeamBlob(blob,beam);
      }
    }
    inline void CleanUp() { return p_rhandler->Reset(); }
  };

}

#endif

