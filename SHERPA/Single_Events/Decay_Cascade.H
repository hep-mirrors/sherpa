#ifndef SHERPA_Single_Events_Decay_Cascade_H
#define SHERPA_Single_Events_Decay_Cascade_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace ATOOLS {
  class Blob_List;
  class Blob;
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

  class Decay_Cascade : public Event_Phase_Handler {
  private :
    const std::vector<PHASIC::Decay_Handler*>& m_decayhandlers;
    SHERPA::Soft_Photon_Handler* p_softphotons;
    int m_mass_smearing, m_qedmode, m_spincorr;
    ATOOLS::Blob_List* p_bloblist;
    
  public :
    Decay_Cascade(const std::vector<PHASIC::Decay_Handler*>& decayhandlers,
                  Soft_Photon_Handler* sph);
    ~Decay_Cascade();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);

    void SetMasses(ATOOLS::Blob* blob, bool usefinalmass);
    void BoostAndStretch(ATOOLS::Blob* blob, const ATOOLS::Vec4D& labmom);

    virtual void TreatInitialBlob(ATOOLS::Blob* blob,
                                  METOOLS::Amplitude2_Tensor* amps);

    METOOLS::Decay_Matrix* FillDecayTree(ATOOLS::Blob * blob, METOOLS::Spin_Density* s0);
    
    bool AttachExtraQED(ATOOLS::Blob* blob,size_t mode=0);
    bool AttachExtraQEDToProductionBlob(ATOOLS::Blob* blob);
    bool AttachExtraQEDRecursively(ATOOLS::Blob* blob, bool aa=false);
    void UpdateDecayBlob(ATOOLS::Blob* blob);
    bool CheckOnshellness(ATOOLS::Blob* blob);

  };
}
#endif

