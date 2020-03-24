#ifndef PHASIC_Channels_Simple_Pole_Channels_H
#define PHASIC_Channels_Simple_Pole_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"

namespace PHASIC {
  class Simple_Pole_RelicDensity: public ISR_Channel_Base {
  protected:
    double m_exponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_sgridkey;
  public:
    Simple_Pole_RelicDensity(const double exponent,const std::string cinfo,
			     ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };
  
  ////////////////////////////////////////////////////////////////////////
  
  class Simple_Pole_Uniform: public ISR_Channel_Base {
  protected:
    double m_exponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Simple_Pole_Uniform(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };
  
  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Forward: public ISR_Channel_Base {
  protected :
    double m_sexponent, m_yexponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Simple_Pole_Forward(const double sexponent,const double yexponent,
			const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Backward: public ISR_Channel_Base {
  protected :
    double m_sexponent, m_yexponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Simple_Pole_Backward(const double sexponent,const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Central: public ISR_Channel_Base {
  protected:
    double m_exponent;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Simple_Pole_Central(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info,const int mode=0);
    void GeneratePoint(const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
  };
}
#endif