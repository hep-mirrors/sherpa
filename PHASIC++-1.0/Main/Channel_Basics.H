#ifndef Channel_Basics_h
#define Channel_Basics_h

#include "Vector.H"

namespace PHASIC {
  class Channel_Basics {
  public:
    /*!
      Rotation of a 4-vector:  p1= rot*p2                          
      According to the flag, either rot (lflag = 0) or p1 (lflag = 1) 
      is calculated.
    */
    static void   Rotat(int lflag,AMATOOLS::Vec4D& p1,AMATOOLS::Vec4D p2,double** rot);
    /*!
      ph is the 4-vector in the rest frame of q                    
      p  is the corresponding 4-vector in the lab frame             
    */
    static void   Boost(int lflag,AMATOOLS::Vec4D q,AMATOOLS::Vec4D& ph,AMATOOLS::Vec4D& p);
    /*!
      This is the squareroot of the ordinary lambda function. 
    */
    static double SqLam(double,double,double);
    /*!
      Selection of ct according to a peaked distribution :
            1/(a + k ct)^cn
      with k = +- 1, and cxm < ct < cxp
    */	    
    static double PeakedDist(double,double,double,double,int,double);
    /*!
      The normalization, i.e. weight, for the peaked distribution.
    */
    static double PeakedWeight(double,double,double,double,int);
    /*!
      Another peaked distribution.
    */
    static double Tj1(double,double,double,double);
    /*!
      And its normalization.
    */
    static double Hj1(double,double,double);
    /*!
      Calculates a minmal pseudo angle for massive particles
      that are treated as massless ones in the phase sapce integration.
      Needed for single t-channelös.
    */
    static double PseudoAngleCut(double,double,double,double); 
  };
}
#endif
