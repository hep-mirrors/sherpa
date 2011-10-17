#ifndef SHERPA_Single_Events_Decay_Handler_Base_H
#define SHERPA_Single_Events_Decay_Handler_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"

namespace ATOOLS {
  class Blob_List;
  class Blob;
  class Cluster_Amplitude;
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

  class Decay_Handler_Base : public ATOOLS::Mass_Selector {
  protected:
    PHASIC::Decay_Map* p_decaymap;
    ATOOLS::Blob_List* p_bloblist;
    ATOOLS::Cluster_Amplitude* p_ampl;
    bool m_extraqed, m_spincorr;
    int m_mass_smearing;

    static bool SortByWidth(ATOOLS::Particle* p1, ATOOLS::Particle* p2);
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
    virtual void AttachExtraQED(ATOOLS::Blob* blob) {}
    virtual void CleanUp();
    ATOOLS::Cluster_Amplitude* ClusterConfiguration(ATOOLS::Blob*const);

    inline void SetBlobList(ATOOLS::Blob_List* bl) { p_bloblist=bl; }
    inline bool SpinCorr() { return m_spincorr; }
    virtual bool Decays(const ATOOLS::Flavour& flav);
    inline PHASIC::Decay_Map* DecayMap() { return p_decaymap; }
  };

}

#endif
