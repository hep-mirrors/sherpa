#ifndef HD_PS_Base_H
#define HD_PS_Base_H

#include "Single_Channel.H"
#include "Multi_Channel.H"
#include "HD_ME_Base.H"  
#include "Flavour.H"
#include "Vector.H"
#include <string>

namespace HADRONS {
  class Hadron_Decay_Channel;

  class HD_Channel_Selector {
  public:
    PHASIC::Single_Channel * GetChannel(int,int,
					const ATOOLS::Flavour *,std::string);
  };
  
  class HD_PS_Base : public PHASIC::Multi_Channel {
  private:
    HD_Channel_Selector  * p_channelselector;
    Hadron_Decay_Channel * p_hdc;
    std::string            m_file;
    double                 m_res,m_error,m_max,m_flux;

    bool     Construct();
    void     AddChannel(std::string,double=1.);
    void     CalculateNormalisedWidth();
    bool     ReadIn();
    bool     WriteOut();
  public:
    HD_PS_Base(Hadron_Decay_Channel *,std::vector<std::string> &,bool &);
    virtual ~HD_PS_Base();
    void     Initialise();
  };
}


#endif
