#ifndef PDF_Main_Cluster_Definitions_Base_H
#define PDF_Main_Cluster_Definitions_Base_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include <vector>

namespace PDF {

  struct CParam {
    double m_kt2, m_op2, m_x, m_mu2;
    int m_kin, m_mode;
    inline CParam(const double &kt2=0.0,
		  const double &op2=0.0,
		  const double &x=0.0,
		  const double &mu2=0.0,
		  const int kin=0,const int mode=0):
      m_kt2(kt2), m_op2(op2), m_x(x), m_mu2(mu2),
      m_kin(kin), m_mode(mode) {}
  };// end of struct CParam

  std::ostream &operator<<(std::ostream &str,const CParam &cp);

  class Cluster_Definitions_Base {
    int m_amode;
  public:

    Cluster_Definitions_Base();
    virtual ~Cluster_Definitions_Base();

    virtual CParam KPerp2
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j=-1,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL,
     const int kin=-1,const int mode=0) = 0;

    virtual ATOOLS::Vec4D_Vector Combine
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL,
     const int kin=0,const int mode=0) = 0;

    virtual CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl);

    int ReCluster(ATOOLS::Cluster_Amplitude *const ampl);

    inline void SetAMode(const size_t mode) { m_amode=mode; }
    inline int AMode() const { return m_amode; }

  };// end of class Cluster_Definitions_Base

}// end of namespace PDF

#endif
