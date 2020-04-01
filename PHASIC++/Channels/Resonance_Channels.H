#ifndef PHASIC_Channels_Resonance_Channels_H
#define PHASIC_Channels_Resonance_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"


namespace PHASIC {
  class Resonance_RelicDensity: public ISR_Channel_Base {
  protected:
    double m_mass, m_width;
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_sgridkey;
  public:
    Resonance_RelicDensity(const double mass,const double width,
			     const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans);
    void GenerateWeight();
  };

  ////////////////////////////////////////////////////////////////////////

  class Resonance_Uniform: public ISR_Channel_Base {
  protected:
    double m_mass, m_width;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Resonance_Uniform(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight();
  };

  ////////////////////////////////////////////////////////////////////////

  class Resonance_Central: public ISR_Channel_Base {
  protected:
    double m_mass, m_width;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Resonance_Central(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight();
  };
  
  ////////////////////////////////////////////////////////////////////////

  class Resonance_Forward: public ISR_Channel_Base {
  protected :
    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Resonance_Forward(const double mass,const double width,const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode);
    void GeneratePoint(const double *rans);
    void GenerateWeight();
  };

  ////////////////////////////////////////////////////////////////////////

  class Resonance_Backward: public ISR_Channel_Base {
  protected :
    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Resonance_Backward(const double mass,const double width,const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info,
		       const size_t mode);
    void GeneratePoint(const double *rans);
    void GenerateWeight();
  };

}
#endif
