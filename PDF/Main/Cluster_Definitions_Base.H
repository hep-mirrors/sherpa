#ifndef PDF_Main_Cluster_Definitions_Base_H
#define PDF_Main_Cluster_Definitions_Base_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include <vector>

namespace PDF {

  struct CParam {
    double m_kt2, m_op2, m_x;
    int m_kin;
    inline CParam(const double &kt2=0.0,
		  const double &op2=0.0,
		  const double &x=0.0,
		  const int kin=0):
      m_kt2(kt2), m_op2(op2), m_x(x), m_kin(kin) {}
  };// end of struct CParam

  std::ostream &operator<<(std::ostream &str,const CParam &cp);

  class Cluster_Definitions_Base {
  public:

    virtual ~Cluster_Definitions_Base();

    virtual CParam KPerp2
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j=-1,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL) = 0;

    virtual ATOOLS::Vec4D_Vector Combine
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL,const int kin=0) = 0;

  };// end of class Cluster_Definitions_Base

}// end of namespace PDF

#endif
