#ifndef Channel_Interface_H
#define Channel_Interface_H

#include "Single_Channel.H"

namespace PHASIC {

  class Channel_Interface : public Single_Channel {
  protected:

    double E,s;
    double E1,E2,Ehat,shat,phat;
    double pt,pl,pt2,pt2min,pt2max,smin,smax;
    double zmin,zmax,z;
    double mass,width;
    int    type;
    int    chnumber;

  public :

    // constructor
    Channel_Interface(int _nin,int _nout,ATOOLS::Flavour *_fl,
		      ATOOLS::Flavour _res=ATOOLS::Flavour(ATOOLS::kf::none)); 

    // setter
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }

    // getters
    void ISRInfo(int &,double &,double &);
    int  ChNumber();

    // member functions
    virtual void GeneratePoint(ATOOLS::Vec4D * p,double *);
    virtual void GenerateWeight(ATOOLS::Vec4D * p);
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *,double * _ran) { GeneratePoint(p,_ran); } 
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *)              { GenerateWeight(p); }

  }; // end of class Channel_Interface

} //end of namespace PHASIC

#endif
