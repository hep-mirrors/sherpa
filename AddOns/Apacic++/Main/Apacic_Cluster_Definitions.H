#ifndef APACIC_Main_Apacic_Cluster_Definitions_H
#define APACIC_Main_Apacic_Cluster_Definitions_H

#include "ATOOLS/Phys/Cluster_Definitions_Base.H"

namespace APACIC {

  class Apacic_Cluster_Definitions: public ATOOLS::Cluster_Definitions_Base {
  public:

    Apacic_Cluster_Definitions();

    ATOOLS::CParam KPerp2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j=-1,int k=-1,
		  const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		  ATOOLS::Mass_Selector *const ms=NULL);

    ATOOLS::Vec4D_Vector Combine
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL);

  };// end of class Apacic_Cluster_Definitions

}// end of namespace APACIC

#endif
