#ifndef REMNANTS_Tools_Kinematics_Generator_H
#define REMNANTS_Tools_Kinematics_Generator_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"

namespace REMNANTS {
  enum class kin_type {
    none   = 0,
    intact = 1,
    coll   = 2,
    DIS1   = 4,
    DIS2   = 5,
    hh     = 8
  };
  class Remnant_Handler;

  class Kinematics_Generator {
  private:
    Remnant_Handler*                  p_rhandler;
    std::array<std::shared_ptr<Remnant_Base>, 2> p_remnants;
    std::array<ATOOLS::Part_List*, 2> p_extracted, p_spectators;
    ATOOLS::Blob*                     p_softblob;

    kin_type m_kintype;

    Primordial_KPerp          m_kperpGenerator;
    ATOOLS::Momenta_Stretcher m_stretcher;

    ATOOLS::Poincare          m_oldcmsboost, m_newcmsboost, m_rotate;
    std::set<ATOOLS::Blob *>  m_boostedblobs;

    double                                m_mass_sum;
    std::array<ATOOLS::Vec4D, 2>          m_inmom, m_checkmom;
    std::array<ATOOLS::ParticleMomMap, 2> m_ktmap;
    ATOOLS::ParticleMomMap                m_shuffledmap;

    unsigned long int m_warns, m_errors;

    void SetKinType();
    bool CollinearKinematics(ATOOLS::Blob_List *blobs);
    bool TransverseKinematics();
    bool TransverseKinematicsHH();
    bool TransverseKinematicsDIS(const size_t & beam);
    bool AdjustShowerInitiators();
    bool BoostConnectedBlob(ATOOLS::Blob * blob,size_t & catcher);
    bool AdjustRemnants();
    bool AdjustFinalStateDIS(const size_t & beam);

    bool CheckHH();
    bool CheckScatter(ATOOLS::Particle * p[2]);
    bool CheckRemnants();
    bool CheckDIS(const size_t & beam);
    const ATOOLS::Vec4D ExtractColourfulFS(
            const size_t& beam, std::vector<ATOOLS::Vec4D>& moms,
            std::vector<double>& masses, std::vector<ATOOLS::Particle*>& parts);
    const ATOOLS::Vec4D ExtractSpectators(
            const size_t& beam, std::vector<ATOOLS::Vec4D>& moms,
            std::vector<double>& masses, std::vector<ATOOLS::Particle*>& parts);

    void InitKTMaps();
    const ATOOLS::Vec4D& ShuffledMomentum(ATOOLS::Particle* part);

  public:
    Kinematics_Generator();
    ~Kinematics_Generator();

    void           Initialize(Remnant_Handler* rhandler);
    ATOOLS::Blob * MakeSoftBlob();
    bool           FillBlobs(ATOOLS::Blob_List *const bloblist);

    inline void SetKinType(const kin_type & type) {
      m_kintype = type;
      if (m_kintype==kin_type::intact) {
	for (size_t beam=0;beam<2;beam++) p_remnants[beam]->SetType(rtp::intact);
      }
    }
    inline Primordial_KPerp * GetKPerp() { return &m_kperpGenerator; }
    void Reset();
  };
}

#endif
