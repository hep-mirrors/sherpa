#ifndef SHERPA_Single_Events_Decay_Handler_Base_H
#define SHERPA_Single_Events_Decay_Handler_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Polarized_CrossSections_Handler.H"

namespace ATOOLS {
  class Blob_List;
  class Blob;
  class Cluster_Amplitude;
  template<typename Scalar>
  class Lorentz_Ten2;
  typedef Lorentz_Ten2<double>  Lorentz_Ten2D;
}

namespace PHASIC {
  class Decay_Map;
  class Decay_Table;
}

namespace METOOLS {
  class Spin_Density;
  class Decay_Matrix;
  class Amplitude2_Tensor;
}

namespace SHERPA {
  class Soft_Photon_Handler;

  class Decay_Handler_Base : public ATOOLS::Mass_Selector {
  protected:
    SHERPA::Soft_Photon_Handler * p_softphotons;
    METOOLS::Polarized_CrossSections_Handler* p_polarization_handler;
    PHASIC::Decay_Map           * p_decaymap;
    ATOOLS::Blob_List           * p_bloblist;
    ATOOLS::Cluster_Amplitude   * p_ampl;
    ATOOLS::Momenta_Stretcher     m_stretcher;
    std::vector<METOOLS::Decay_Matrix> m_decaymatrices;
    bool m_spincorr, m_decaychainend, m_cluster, m_polcrosssec;
    int m_mass_smearing;
    size_t m_qedmode;
    unsigned long int m_oserrors;
    bool m_specialtauspincorr;
    bool DoSpecialDecayTauSC(ATOOLS::Particle* part);
  public:
    Decay_Handler_Base();
    ~Decay_Handler_Base();
    void SetMasses(ATOOLS::Blob* blob, bool usefinalmass);
    void BoostAndStretch(ATOOLS::Blob* blob, const ATOOLS::Vec4D& labmom);

    // have to reimplement for hadrons, lund:
    virtual void TreatInitialBlob(ATOOLS::Blob* blob,
                                  METOOLS::Amplitude2_Tensor* amps,
                                  const ATOOLS::Particle_Vector& origparts);
    // have to reimplement for lund:
    virtual bool DiceMass(ATOOLS::Particle* part, double max);
    virtual METOOLS::Amplitude2_Tensor*
    FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);
    virtual METOOLS::Decay_Matrix*
    FillDecayTree(ATOOLS::Blob * blob, METOOLS::Spin_Density* s0);
    virtual void CreateDecayBlob(ATOOLS::Particle* inpart)=0;
    virtual void UpdateDecayBlob(ATOOLS::Blob* blob);
    virtual bool CheckOnshellness(ATOOLS::Blob* blob);
    virtual void CleanUp();
    ATOOLS::Cluster_Amplitude* ClusterConfiguration(ATOOLS::Blob*const);

    inline void SetBlobList(ATOOLS::Blob_List* bl) { p_bloblist=bl; }
    inline std::vector<METOOLS::Decay_Matrix> GetDecayMatrices() { return m_decaymatrices; }
    inline bool SpinCorr() { return m_spincorr; }
    inline bool PolCrossSec(){ return  m_polcrosssec; }
    inline METOOLS::Polarized_CrossSections_Handler* GetPolarizationHandler() { return p_polarization_handler; }
    virtual bool Decays(const ATOOLS::Flavour& flav);
    virtual bool CanDecay(const ATOOLS::Flavour& flav);
    inline PHASIC::Decay_Map* DecayMap() { return p_decaymap; }

    inline void SetSoftPhotonHandler(Soft_Photon_Handler* sph)
    {  p_softphotons=sph; }
    inline Soft_Photon_Handler* GetSoftPhotonHandler()
    { return p_softphotons; }
  };


}

#endif
