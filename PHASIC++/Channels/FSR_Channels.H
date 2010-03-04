#ifndef PHASIC_Channels_FSR_Channels_H
#define PHASIC_Channels_FSR_Channels_H

#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {

  class Phase_Space_Handler;

  class FSR_Channels: public Multi_Channel {
  private:

    Phase_Space_Handler *p_psh;

    int m_inttype;

    void DropRedundantChannels();
    bool CompareCh(std::string C1,std::string C2);
    bool Compare(const ATOOLS::Vec4D *p1,const ATOOLS::Vec4D *p2);

  public:

    inline FSR_Channels(Phase_Space_Handler *const psh,
			const std::string &name):
      Multi_Channel(name), p_psh(psh) {}

    bool Initialize();

  };// end of class FSR_Channels

}//end of namespace PHASIC

#endif
