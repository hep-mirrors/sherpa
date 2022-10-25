#ifndef SHERPA_SoftPhysics_Soft_Collision_Handler_H
#define SHERPA_SoftPhysics_Soft_Collision_Handler_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "PDF/Main/Cluster_Definitions_Base.H"

namespace ATOOLS {
  class Cluster_Amplitude;
  class Cluster_Definitions_Base;
}
namespace SHRIMPS { class Shrimps; }
namespace AMISIC  { class Amisic; }
namespace PDF     { class ISR_Handler; }
namespace BEAM    { class Beam_Spectra_Handler; }
namespace MODEL   { class Model_Base;           }

namespace SHERPA {
  class Soft_Collision_Handler: public ATOOLS::Terminator_Object {
  private:
    struct scmode {
      enum code {
	none,
	shrimps,
	amisic
      };
    };
    std::string  m_scmodel;
    std::string  m_dir,m_sfile;
    scmode::code m_mode;

    SHRIMPS::Shrimps           * p_shrimps;
    AMISIC::Amisic             * p_amisic;
    void PrepareTerminate(); 
  public:
    Soft_Collision_Handler(AMISIC::Amisic * amisic,
			   SHRIMPS::Shrimps * shrimps);
    ~Soft_Collision_Handler();

    ATOOLS::Return_Value::code  GenerateMinimumBiasEvent(ATOOLS::Blob_List * blobs);
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob *const bl);
    
    void SetShrimps(SHRIMPS::Shrimps * shrimps=NULL);
    void CleanUp();
    
    inline void SetAmisic(AMISIC::Amisic * amisic) {
      p_amisic  = amisic;
      m_scmodel = std::string("Amisic");
      m_mode    = scmode::amisic;
    }    
    inline const std::string &  Soft_CollisionModel() const { 
      return m_scmodel; 
    }
    inline SHRIMPS::Shrimps * GetShrimps() const { 
      return p_shrimps; 
    }
    inline AMISIC::Amisic * GetAmisic() const { 
      return p_amisic; 
    }
  };// end of class Soft_Collision_Handler
}
#endif

