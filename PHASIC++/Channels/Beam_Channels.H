#ifndef PHASIC_Channels_Beam_Channels_H
#define PHASIC_Channels_Beam_Channels_H

#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {
  class Phase_Space_Handler;

  class Beam_Channels: public Multi_Channel {
  private:
    BEAM::beammode::code m_beammode;
    channelinfos         m_beamparams;

    
    Phase_Space_Handler *p_psh;
    std::string          m_keyid;
    
    bool MakeChannels();
    bool CreateChannels();
    void AddSimplePole(const size_t & beamno,const int & beam);
    void AddResonance(const size_t & beamno,const int & beam);
    void AddThreshold(const size_t & beamno,const int & beam);
    void AddLaserBackscattering(const size_t & beamno,const int & beam);    
  public:
    Beam_Channels(Phase_Space_Handler *const psh,const std::string &name);
    ~Beam_Channels() {}
    bool Initialize();
    inline channelinfos &Parameters()  { return m_beamparams;  } 
  };// end of class Beam_Channels

}//end of namespace PHASIC

#endif
