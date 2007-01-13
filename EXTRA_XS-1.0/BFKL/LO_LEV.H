#ifndef LO_LEV_H
#define LO_LEV_H

#include "LEV_Base.H"

namespace EXTRAXS {

  class G_GG: public LEV_Base {
  public:
    
    // constructor
    G_GG();

    // member functions
    double Value(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;
    double MajorValue(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		      const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;

    double MajorIntegral();

  };// end of class LEV_Base

}// end of namespace EXTRAXS

#endif
