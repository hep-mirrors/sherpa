#ifndef Channel_Interface_H
#define Channel_Interface_H

#include "Single_Channel.H"

namespace PHASIC {

  class Channel_Interface : public Single_Channel {
  protected:

    double E, Ehat, s, smax, pt2, pt2max, pt2min, pt, mass, width;

    int type, chnumber;

  public :

    // constructor
    Channel_Interface(int nin,int nout,ATOOLS::Flavour *flavour,
		      ATOOLS::Flavour res=ATOOLS::Flavour(ATOOLS::kf::none)); 

    // setter
    void SetChNumber(int _chnumber) { chnumber = _chnumber; }

    // getters
    void ISRInfo(int &type,double &mass,double &width);
    int  ChNumber();

    // member functions
    virtual void GeneratePoint(ATOOLS::Vec4D *p,double *ran);
    virtual void GenerateWeight(ATOOLS::Vec4D *p);

    inline void GeneratePoint(ATOOLS::Vec4D *p,
			      ATOOLS::Cut_Data *cuts,double *ran)       { GeneratePoint(p,ran); } 
    inline void GenerateWeight(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts) { GenerateWeight(p);    }

  };// end of class Channel_Interface

}//end of namespace PHASIC

#endif
