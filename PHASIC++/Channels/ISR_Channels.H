#ifndef PHASIC_Channels_ISR_Channels_H
#define PHASIC_Channels_ISR_Channels_H

#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {

  class Phase_Space_Handler;

  class ISR_Channels: public Multi_Channel {
  private:

    channelinfos m_isrparams;

    Phase_Space_Handler *p_psh;

    bool MakeChannels();
    bool CreateChannels();

  public:

    inline ISR_Channels(Phase_Space_Handler *const psh,
			const std::string &name):
      Multi_Channel(name), p_psh(psh) {}

    bool Initialize();

    inline channelinfos &Parameters()  { return m_isrparams;  } 

  };// end of class ISR_Channels

}//end of namespace PHASIC

#endif
