#ifndef SHERPA_SoftPhysics_Soft_Collision_Handler_H
#define SHERPA_SoftPhysics_Soft_Collision_Handler_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"

namespace ATOOLS {
  class Cluster_Amplitude;
  class Cluster_Definitions_Base;
}
namespace SHRIMPS { class Shrimps; }
namespace AMISIC  { class Amisic; }
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
    bool         m_bunch_rescatter;
    std::string  m_scmodel;
    std::string  m_dir,m_sfile;
    scmode::code m_mode;

    SHRIMPS::Shrimps           * p_shrimps;
    AMISIC::Amisic             * p_amisic;

    void   PrepareTerminate();
    double ImpactParameter();
  public:
    Soft_Collision_Handler(AMISIC::Amisic * amisic,SHRIMPS::Shrimps * shrimps,
			   const bool bunch_rescatter=false);
    ~Soft_Collision_Handler();

    ATOOLS::Return_Value::code  GenerateMinimumBiasEvent(ATOOLS::Blob_List * blobs);
    void                        SetPosition(const size_t & beam,const ATOOLS::Vec4D & pos);
    ATOOLS::Return_Value::code  GenerateBunchRescatter(ATOOLS::Blob_List * blobs);
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
    inline const bool & BeamRescatter() { return m_bunch_rescatter; }
  };// end of class Soft_Collision_Handler
  typedef std::map<PDF::isr::id, Soft_Collision_Handler *> Soft_Collision_Handler_Map;
}
#endif

