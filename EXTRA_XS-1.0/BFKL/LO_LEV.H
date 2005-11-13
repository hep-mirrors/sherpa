#ifndef LO_LEV_H
#define LO_LEV_H

#include "LEV_Base.H"

namespace EXTRAXS {

  class G_GG: public LEV_Base {
  public:
    
    // constructor
    G_GG();

    // member functions
    double Value(const double &y,const double &kt) const;
    double MajorValue(const double &y,const double &kt) const;

    double MajorIntegral();

  };// end of class LEV_Base

}// end of namespace EXTRAXS

#endif
