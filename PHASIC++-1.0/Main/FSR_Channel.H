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
    double mass,width;
    int    type;
    int    chnumber;
  public :
    S1Channel(int,int,APHYTOOLS::Flavour *,
	      APHYTOOLS::Flavour=APHYTOOLS::Flavour(APHYTOOLS::kf::none)); 
    void GeneratePoint(AMATOOLS::Vec4D * p,double *);
    void GenerateWeight(AMATOOLS::Vec4D * p);
    void GeneratePoint(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *)              { GenerateWeight(p); }
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };


  class T1Channel : public Single_Channel {
    /*!
      This channel is meant to be the default t-channel integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
    double mass,width;
    int    type;
    int    chnumber;
  public :
    T1Channel(int,int,APHYTOOLS::Flavour *,
	      APHYTOOLS::Flavour=APHYTOOLS::Flavour(APHYTOOLS::kf::none)); 
    void GeneratePoint(AMATOOLS::Vec4D * p,double *);
    void GenerateWeight(AMATOOLS::Vec4D * p);
    void GeneratePoint(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *)              { GenerateWeight(p); }
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };



  class U1Channel : public Single_Channel {
    /*!
      This channel is meant to be the default u-channel integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
    double mass,width;
    int    type;
    int    chnumber;
  public :
    U1Channel(int,int,APHYTOOLS::Flavour *,
	      APHYTOOLS::Flavour=APHYTOOLS::Flavour(APHYTOOLS::kf::none)); 
    void GeneratePoint(AMATOOLS::Vec4D * p,double *);
    void GenerateWeight(AMATOOLS::Vec4D * p);
    void GeneratePoint(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(AMATOOLS::Vec4D * p,APHYTOOLS::Cut_Data *)              { GenerateWeight(p); }
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };
}


#endif
