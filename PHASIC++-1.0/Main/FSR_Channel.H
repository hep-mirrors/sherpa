#ifndef FSR_Channel_h
#define FSR_Channel_h

#include "Single_Channel.H"

namespace PHASIC {
  class S1Channel : public Single_Channel {
    /*!
      This channel is meant to be the default s-channel integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
  public :
    S1Channel(int,int,APHYTOOLS::Flavour *); 
    void GeneratePoint(AMATOOLS::vec4d * p,double *);
    void GenerateWeight(AMATOOLS::vec4d * p);
  };


  class T1Channel : public Single_Channel {
    /*!
      This channel is meant to be the default s-channel integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
  public :
    T1Channel(int,int,APHYTOOLS::Flavour *); 
    void GeneratePoint(AMATOOLS::vec4d * p,double *);
    void GenerateWeight(AMATOOLS::vec4d * p);
  };



  class U1Channel : public Single_Channel {
    /*!
      This channel is meant to be the default s-channel integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
  public :
    U1Channel(int,int,APHYTOOLS::Flavour *); 
    void GeneratePoint(AMATOOLS::vec4d * p,double *);
    void GenerateWeight(AMATOOLS::vec4d * p);
  };
}


#endif
