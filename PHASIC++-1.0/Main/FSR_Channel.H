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
    S1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(ATOOLS::kf::none)); 
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
  };

  class SimpleQCDChannel : public Single_Channel {
  public :
    SimpleQCDChannel(int nin,int nout,ATOOLS::Flavour *flavs); 
    void GeneratePoint(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts,double *ran);
    void GenerateWeight(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts);
    void ISRInfo(int &,double &,double &);
    std::string ChID();
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
    T1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(ATOOLS::kf::none)); 
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
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
    U1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(ATOOLS::kf::none)); 
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
  };


  class Decay2Channel : public Single_Channel {
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
    Decay2Channel(int,int,ATOOLS::Flavour *,
		  ATOOLS::Flavour=ATOOLS::Flavour(ATOOLS::kf::none)); 
    void GeneratePoint(ATOOLS::Vec4D * p,double *);
    void GenerateWeight(ATOOLS::Vec4D * p);
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *)              { GenerateWeight(p); }
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };
}


#endif
