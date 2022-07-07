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
namespace SHRIMPS {
  class Shrimps;
  class Cluster_Algorithm;
}
namespace AMISIC {
  class Amisic;
}
namespace PDF    { class ISR_Handler; }
namespace BEAM   { class Beam_Spectra_Handler; }
namespace MODEL  { class Model_Base;           }

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
    std::string  m_dir,m_sfile;
    std::string  m_softcollisionmodel;
    scmode::code m_mode;

    SHRIMPS::Shrimps           * p_shrimps;
    SHRIMPS::Cluster_Algorithm * p_cluster;
    AMISIC::Amisic             * p_amisic;
    void PrepareTerminate(); 
  public:
    Soft_Collision_Handler(MODEL::Model_Base *model,
			   BEAM::Beam_Spectra_Handler*,
                           PDF::ISR_Handler*);
    ~Soft_Collision_Handler();

    ATOOLS::Return_Value::code
    GenerateMinimumBiasEvent(ATOOLS::Blob_List* blobs);
    void CleanUp();
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob *const bl);

    inline void SetAmisic(AMISIC::Amisic * amisic) {
      p_amisic             = amisic;
      m_softcollisionmodel = std::string("Amisic");
      m_mode               = scmode::amisic;
    }
    inline void 
    SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {
      //p_cluster->SetClusterDefinitions(defs);
    }
    
    inline const std::string &  Soft_CollisionModel() const { 
      return m_softcollisionmodel; 
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

