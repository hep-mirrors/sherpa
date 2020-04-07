#ifndef HADRONS_Main_Hadron_Decay_Handler_H
#define HADRONS_Main_Hadron_Decay_Handler_H

#include <set>
#include <iterator>

#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "PHASIC++/Decays/Decay_Handler.H"


namespace ATOOLS { 
  class Blob;
  class Blob_List;
}

namespace HADRONS {
  class Mixing_Handler;

  class Hadron_Decay_Handler : public PHASIC::Decay_Handler {
    HADRONS::Mixing_Handler* p_mixinghandler;

    long int bs, bs_cs_semilep, bs_ccss, bs_ccds, bs_cuds;
    long int bs_cs_semilep_rej, bs_ccss_rej, bs_ccds_rej, bs_cuds_rej;

    bool RejectExclusiveChannelsFromFragmentation(ATOOLS::Blob*);
    void MakeShowerBlobsForPartonicDecays(ATOOLS::Blob* initialblob, ATOOLS::Blob_List* bloblist);

    void SetPosition(ATOOLS::Blob* blob);

    void RegisterSettings();
  public :
    Hadron_Decay_Handler();
    ~Hadron_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Blob_List* bloblist, ATOOLS::Particle* part);
    bool VetoDecayAndPrepForNew(ATOOLS::Blob* blob);
    void AfterTreatInitialBlob(ATOOLS::Blob* blob, ATOOLS::Blob_List* bloblist);
    void BeforeFillDecayTree(ATOOLS::Blob * blob);
    METOOLS::Amplitude2_Tensor* FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);

    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }
    inline std::string Name() { return "HADRONS++"; }
  };

}

#endif
