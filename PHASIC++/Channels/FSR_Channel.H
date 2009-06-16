#ifndef PHASIC_Channels_FSR_Channel_H
#define PHASIC_Channels_FSR_Channel_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Vegas.H"

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
    Vegas *p_vegas;
  public :
    S1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~S1Channel();
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
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
    Vegas *p_vegas;
  public :
    T1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~T1Channel();
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
  };


  class T2Channel : public Single_Channel {
    /*!
      This channel is for bfkl stuff only.
      It dices the momenta in the ladder 
      from bottom to top.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
    double mass,width;
    int    type;
    int    chnumber;
    Vegas *p_vegas;
  public :
    T2Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~T2Channel();
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
  };

  class T3Channel : public Single_Channel {
    /*!
      This channel is for bfkl stuff only.
      It dices the momenta in the ladder 
      from top to bottom.
    */
    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
    double mass,width;
    int    type;
    int    chnumber;
    Vegas *p_vegas;
  public :
    T3Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~T3Channel();
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
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
    Vegas *p_vegas;
  public :
    U1Channel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~U1Channel();
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
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
    Decay2Channel(int,int,const ATOOLS::Flavour *,
		  ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    void GeneratePoint(ATOOLS::Vec4D * p,double *);
    void GenerateWeight(ATOOLS::Vec4D * p);
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *)              { GenerateWeight(p); }
    void ISRInfo(int &,double &,double &);
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };
}


#endif
