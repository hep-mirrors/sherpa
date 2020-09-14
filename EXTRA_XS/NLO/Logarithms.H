#ifndef EXTRAXS__NLO__Logarithms_H
#define EXTRAXS__NLO__Logarithms_H

#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;

namespace EXTRAXS {

  inline double HTheta(const double &x)
  {
    return x>0.0?1.0:0.0;
  }
  
  inline Complex LnRat(const double &x,const double &y)
  {
    return log(dabs(x/y))-Complex(0.0,M_PI)*(HTheta(-x)-HTheta(-y));
  }

  Complex L0(const double &x,const double &y);
  Complex L1(const double &x,const double &y);
  Complex L2(const double &x,const double &y);

  Complex Ls0(const double &x1,const double &y1,
	      const double &x2,const double &y2);
  Complex Ls1(const double &x1,const double &y1,
	      const double &x2,const double &y2);
  Complex Lsm1(const double &x1,const double &y1,
	       const double &x2,const double &y2);

}// end of namespace EXTRAXS

#endif
