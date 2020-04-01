#ifndef PHASIC_Channels_ISR_Channels_H
#define PHASIC_Channels_ISR_Channels_H

#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {
  class Phase_Space_Handler;
  
  class ISR_Channels: public Multi_Channel {
  private:
    Phase_Space_Handler * p_psh;
    std::string           m_keyid;
    PDF::ISR_Handler    * p_isrhandler;
    PDF::isrmode::code    m_isrmode;
    PDF::isrtype::code    m_isrtype[2];
    channelinfos          m_isrparams;
    std::set<double>      m_yexponents, m_spexponents;

    bool MakeChannels();
    void CheckForStructuresFromME();
    bool CreateChannels();
    void AddSimplePole(const size_t & chno,const size_t & mode);
    void AddResonance(const size_t & chno,const size_t & mode);
    void AddThreshold(const size_t & chno,const size_t & mode);
    void AddLeadingLog(const size_t & chno,const size_t & mode);    
  public:
    ISR_Channels(Phase_Space_Handler *const psh,const std::string &name);
    inline ~ISR_Channels() {}
    inline bool Initialize()          { return MakeChannels(); }
    inline channelinfos &Parameters() { return m_isrparams;  } 
  };// end of class ISR_Channels

}//end of namespace PHASIC

#endif
