#ifndef PHASIC_Channels_Beam_Channels_H
#define PHASIC_Channels_Beam_Channels_H

#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"

namespace PHASIC {
  class Phase_Space_Handler;

  class Beam_Channels: public Multi_Channel {
  private:
    Phase_Space_Handler        * p_psh;
    std::string                  m_keyid;
    BEAM::Beam_Spectra_Handler * p_beamspectra;
    BEAM::beammode               m_beammode;
    BEAM::beamspectrum           m_beamtype[2];
    channelinfos                 m_beamparams;
    std::set<double>             m_yexponents;


    bool MakeChannels();
    bool DefineColliderChannels();
    void CheckForStructuresFromME();
    bool CreateChannels();
    void AddSimplePole(const size_t & chno,const size_t & mode);
    void AddResonance(const size_t & chno,const size_t & mode);
    void AddThreshold(const size_t & chno,const size_t & mode);
    void AddLaserBackscattering(const size_t & chno,const size_t & mode);
    void AddExponential(const size_t & chno,const size_t & mode);
  public:
    Beam_Channels(Phase_Space_Handler *const psh,const std::string &name);
    ~Beam_Channels() {}
    bool Initialize();
    inline channelinfos &Parameters()  { return m_beamparams;  }
  };// end of class Beam_Channels

}//end of namespace PHASIC

#endif
