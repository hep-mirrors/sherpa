#ifndef PHASIC_Channels_LBS_Compton_Peak_Channels_H
#define PHASIC_Channels_LBS_Compton_Peak_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"


namespace PHASIC {
  class LBS_Compton_Peak_Uniform: public ISR_Channel_Base {
  protected:
    double m_exponent, m_pole;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    LBS_Compton_Peak_Uniform(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info,
			     const int mode=0);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };
  
  class LBS_Compton_Peak_Forward: public ISR_Channel_Base {
  protected :
    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    LBS_Compton_Peak_Forward(const double exponent,const double pole,const double yexponent,
			     const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };

  class LBS_Compton_Peak_Backward: public ISR_Channel_Base {
  protected :
    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    LBS_Compton_Peak_Backward(const double exponent,const double pole,const double yexponent,
			      const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };

  class LBS_Compton_Peak_Central: public ISR_Channel_Base {
  protected:
    double m_exponent, m_pole;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    LBS_Compton_Peak_Central(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info,
			     const int mode=0);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };
}
#endif
