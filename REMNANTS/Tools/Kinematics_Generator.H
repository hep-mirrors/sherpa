#ifndef REMNANTS_Tools_Kinematics_Generator_H
#define REMNANTS_Tools_Kinematics_Generator_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"

namespace REMNANTS {
  struct kin_type {
    enum code {
      none   = 0,
      intact = 1,
      coll   = 2,
      DIS1   = 4,
      DIS2   = 5,
      hh     = 8
    };
  };
  class Remnant_Handler;
  
  class Kinematics_Generator {
  private:
    Remnant_Handler   * p_rhandler;
    Remnant_Base      * p_remnants[2];
    ATOOLS::Part_List * p_extracted[2], * p_spectators[2];
    ATOOLS::Blob      * p_softblob;
    
    kin_type::code  m_kintype;
    
    Primordial_KPerp          m_kperpGenerator; 
    ATOOLS::Momenta_Stretcher m_stretcher;

    ATOOLS::Poincare          m_oldcmsboost, m_newcmsboost, m_rotate;
    std::set<ATOOLS::Blob *>  m_boostedblobs;

    ATOOLS::Vec4D  m_inmom[2], m_checkmom[2];
    ParticleMomMap m_ktmap, m_shuffledmap;

    void SetKinType(Remnant_Handler * const rhandler);
    bool CollinearKinematics();
    bool TransverseKinematics();
    void InitKTMaps();
    bool AdjustShowerInitiators();
    bool ConstructMomentaAndBoosts(ATOOLS::Particle * part[2]);
    bool BoostConnectedBlob(ATOOLS::Blob * blob,size_t & catcher);
    bool AdjustRemnants();

    bool Check();
    bool CheckScatter(ATOOLS::Particle * p[2]);
    bool CheckRemnants();
    
    const ATOOLS::Vec4D & ShuffledMomentum(ATOOLS::Particle *const part);
  public:
    Kinematics_Generator();
    ~Kinematics_Generator();

    void Initialize(Remnant_Handler * const rhandler,
		    const std::string & path, const std::string & file);
    ATOOLS::Blob * MakeSoftBlob();
    bool FillBlobs(ATOOLS::Blob_List *const bloblist);
      
    void Reset();
  };
}

#endif
