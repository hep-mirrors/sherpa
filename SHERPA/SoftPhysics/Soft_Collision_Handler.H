#ifndef SHERPA_SoftPhysics_Soft_Collision_Handler_H
#define SHERPA_SoftPhysics_Soft_Collision_Handler_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Exception.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "SHERPA/SoftPhysics/Cluster_Algorithm.H"

namespace ATOOLS { class Cluster_Amplitude; class Cluster_Definitions_Base; }
namespace PDF    { class ISR_Handler; }
namespace BEAM   { class Beam_Spectra_Handler; }
namespace SHRIMPS  { class Shrimps; }

namespace SHERPA {
  class Soft_Collision_Handler: public ATOOLS::Terminator_Object {
  private:
    std::string m_dir,m_file,m_sfile;
    std::string m_softcollisionmodel;
    int         m_mode;

    SHRIMPS::Shrimps      * p_shrimps;
    Cluster_Algorithm   m_cluster;
    void PrepareTerminate(); 
  public:
    Soft_Collision_Handler(std::string path,std::string file,
			   BEAM::Beam_Spectra_Handler *beam,
			   PDF::ISR_Handler *isr);
    ~Soft_Collision_Handler();

    ATOOLS::Return_Value::code 
    GenerateMinimumBiasEvent(ATOOLS::Blob_List * blobs,
			     double & weight);
    void CleanUp();
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob *const bl);
    
    inline void 
    SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {
      m_cluster.SetClusterDefinitions(defs);
    }
    inline const std::string &  Soft_CollisionModel() const { 
      return m_softcollisionmodel; 
    }
    inline SHRIMPS::Shrimps * GetShrimps() const { 
      return p_shrimps; 
    }
  };// end of class Soft_Collision_Handler
}
#endif

