#ifndef PHASIC__Scales__Cluster_Definitions_H
#define PHASIC__Scales__Cluster_Definitions_H

#include "PDF/Main/Cluster_Definitions_Base.H"

namespace PHASIC {

  class Cluster_Definitions: public PDF::Cluster_Definitions_Base {
  private:

    int m_kfac, m_nproc, m_pstype;

  public:

    inline Cluster_Definitions
    (const int kfac,const int nproc,const int pstype):
      m_kfac(kfac), m_nproc(nproc), m_pstype(pstype) {}

    PDF::Cluster_Param Cluster(const PDF::Cluster_Config &ca);

  };// end of class Cluster_Definitions

}// end of namespace PHASIC

#endif
