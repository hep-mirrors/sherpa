#ifndef ATOOLS_Phys_Cluster_Definitions_Base_H
#define ATOOLS_Phys_Cluster_Definitions_Base_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include <vector>

namespace ATOOLS {

  struct CParam {
    double m_kt2, m_op2, m_x;
    inline CParam(const double &kt2=0.0,
		  const double &op2=0.0,
		  const double &x=0.0):
      m_kt2(kt2), m_op2(op2), m_x(x) {}
  };// end of struct CParam

  std::ostream &operator<<(std::ostream &str,const CParam &cp);

  class Cluster_Definitions_Base {
  public:

    virtual ~Cluster_Definitions_Base();

    virtual CParam KPerp2(const Cluster_Amplitude &ampl,
			  int i,int j=-1,int k=-1,
			  const Flavour &mo=Flavour(kf_none),
			  Mass_Selector *const ms=NULL) = 0;

    virtual Vec4D_Vector Combine(const Cluster_Amplitude &ampl,
				 int i,int j,int k=-1,
				 const Flavour &mo=Flavour(kf_none),
				 Mass_Selector *const ms=NULL) = 0;

  };// end of class Cluster_Definitions_Base

}// end of namespace ATOOLS

#endif
