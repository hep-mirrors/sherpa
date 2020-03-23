#ifndef PHASIC_Channels_ISR_Channels_H
#define PHASIC_Channels_ISR_Channels_H

#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {
  class Phase_Space_Handler;
  
  class ISR_Channels: public Multi_Channel {
  private:
    Phase_Space_Handler * p_psh;
    channelinfos          m_isrparams;
    std::string           m_keyid;
    
    bool MakeChannels();
    bool CreateChannels();
  public:
    inline ISR_Channels(Phase_Space_Handler *const psh,const std::string &name) :
      p_psh(psh), Multi_Channel(name), m_keyid("ISR") {}
    inline ~ISR_Channels() {}
    inline bool Initialize()          { return MakeChannels(); }
    inline channelinfos &Parameters() { return m_isrparams;  } 
  };// end of class ISR_Channels

}//end of namespace PHASIC

#endif
