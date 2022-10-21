#ifndef PHASIC_Channels_ISR_Channel_Interface_H
#define PHASIC_Channels_ISR_Channel_Interface_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"

namespace PHASIC {

  class Multi_Channel;
  
  class ISR_Channel_Interface: public ISR_Channel_Base {
  protected:

    Multi_Channel *p_fsr;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;

  public:

    ISR_Channel_Interface(Multi_Channel *const fsr,
			  const std::string &cinfo,
			  ATOOLS::Integration_Info *info);

    void GeneratePoint(const double *rans);
    void GenerateWeight(const int &mode=0);
    void AddPoint(double);

  };// end of class ISR_Channel_Interface

}

#endif

    
