#ifndef Channel_Basics_H
#define Channel_Basics_H

#include "Vector.H"

namespace PHASIC {

  class Channel_Basics {
  public:

    static void   Rotat(int lflag,ATOOLS::Vec4D& p1,ATOOLS::Vec4D p2,double** rot);
    static void   Boost(int lflag,ATOOLS::Vec4D q,ATOOLS::Vec4D& ph,ATOOLS::Vec4D& p);
    static double SqLam(double,double,double);
    static double PeakedDist(double,double,double,double,int,double);
    static double PeakedWeight(double,double,double,double,int);
    static double Tj1(double,double,double,double);
    static double Hj1(double,double,double);
    static double PseudoAngleCut(double,double,double,double); 

  };// end of class Channel_Basics

  /*!
    Rotation of a 4-vector:  p1= rot*p2                          
    According to the flag, either rot (lflag = 0) or p1 (lflag = 1) 
    is calculated.
  */

  /*!
    ph is the 4-vector in the rest frame of q                    
    p  is the corresponding 4-vector in the lab frame             
  */

  /*!
    This is the squareroot of the ordinary lambda function. 
  */

  /*!
    Selection of ct according to a peaked distribution :
    1/(a + k ct)^cn
    with k = +- 1, and cxm < ct < cxp
  */	    

  /*!
    The normalization, i.e. weight, for the peaked distribution.
  */

  /*!
    Another peaked distribution.
  */

  /*!
    And its normalization.
  */

  /*!
    Calculates a minmal pseudo angle for massive particles
    that are treated as massless ones in the phase sapce integration.
    Needed for single t-channelös.
  */

}// end of namespace PHASIC

#endif
