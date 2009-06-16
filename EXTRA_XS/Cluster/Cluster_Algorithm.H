#ifndef EXTRA_XS_Cluster_Cluster_Algorithm_H
#define EXTRA_XS_Cluster_Cluster_Algorithm_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace ATOOLS { class Cluster_Definitions_Base; }

namespace EXTRAXS {

  class Single_Process;

  class Cluster_Algorithm {
  private:

    ATOOLS::Cluster_Amplitude        *p_ampl;
    ATOOLS::Cluster_Definitions_Base *p_clus;

    bool SetColors();

  public:

    Cluster_Algorithm();

    ~Cluster_Algorithm();

    bool Cluster(Single_Process *const xs);
    
    // inline functions
    inline ATOOLS::Cluster_Amplitude *Amplitude() { return p_ampl; }

    inline void SetClusterDefinitions
    (ATOOLS::Cluster_Definitions_Base *const cb) { p_clus=cb; }

  };// end of class Cluster_Algorithm

}// end of namespace EXTRAXS

#endif
