#ifndef Two_Body_PSs_h
#define Two_Body_PSs_h

#include "Single_Channel.H"
#include "Vector.H"

namespace HADRONS {
  class Iso2Channel : public PHASIC::Single_Channel {
    /*!
      This channel is meant to be the default 2-body integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    ATOOLS::Vec4D m_decvec;
    int           m_chnumber;
  public :
    Iso2Channel(const ATOOLS::Flavour *);
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL,double * _ran=NULL);
    void GeneratePoint(ATOOLS::Vec4D * p,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL);
    void GenerateWeight(ATOOLS::Vec4D * p);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return std::string("Iso2");}
  };
}

#endif
