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
    void MPISync() { p_vegas->MPISync(); }
    void Optimize() { p_vegas->Optimize(); }
    void EndOptimize() { p_vegas->EndOptimize(); }
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId) { p_vegas->ReadIn(pId); }
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
    void MPISync() { p_vegas->MPISync(); }
    void Optimize() { p_vegas->Optimize(); }
    void EndOptimize() { p_vegas->EndOptimize(); }
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId) { p_vegas->ReadIn(pId); }
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
    void MPISync() { p_vegas->MPISync(); }
    void Optimize() { p_vegas->Optimize(); }
    void EndOptimize() { p_vegas->EndOptimize(); }
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId) { p_vegas->ReadIn(pId); }
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
    void MPISync() { }
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
  };


  class NoChannel : public Single_Channel {
    /*!
      This channel is meant to be the integrator for s-channel instanton production
      and similar processes, where an s-channel object with no specific mass is
produced.
    */
    double E,s;
    double pt2min,pt2max,smin,smax;
    double mass,width;
    int    type;
    int    chnumber;
  public :
    NoChannel(int,int,ATOOLS::Flavour *,
	      ATOOLS::Flavour=ATOOLS::Flavour(kf_none)); 
    ~NoChannel() {}
    void GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *,double * _ran);
    void GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *);
    void ISRInfo(int &,double &,double &);
    void MPISync()                  { }
    void Optimize()                 { }
    void EndOptimize()              { }
    void WriteOut(std::string pId)  { }
    void ReadIn(std::string pId)    { }
    int  ChNumber()                 { return chnumber;      }
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }
    std::string ChID();
  };
}


#endif
