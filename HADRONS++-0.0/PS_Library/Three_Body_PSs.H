#ifndef Three_Body_PSs_h
#define Three_Body_PSs_h

#include "Single_Channel.H"
#include "Vector.H"

namespace HADRONS {
  class Dalitz : public PHASIC::Single_Channel {
    /*!
      This channel is meant to be the default 2-body integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.
    */
    ATOOLS::Vec4D m_decvec, m_pvec;
    double        m_pmass, m_pwidth, m_sexp, m_smin, m_smax;
    int           m_p1, m_p2, m_dir, m_mode;
    std::string   m_propmass;
    int           m_chnumber;
  public :
    Dalitz(const ATOOLS::Flavour *,const ATOOLS::Flavour,const int,const int);
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL,double * _ran=NULL);
    void GeneratePoint(ATOOLS::Vec4D * p,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL);
    void GenerateWeight(ATOOLS::Vec4D * p);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return name;}
  };
}

#endif
