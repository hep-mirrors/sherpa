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
    void GeneratePoint(const double *rans);
      void GenerateWeight(const int & mode=0);
  };

  class Simple_Pole_DM_Annihilation: public ISR_Channel_Base {
  protected:
    double m_exponent, m_mass[2];
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_sgridkey, m_xkey, m_xgridkey, m_cosxikey, m_cosgridkey;
  public:
    Simple_Pole_DM_Annihilation(const double exponent, const double mass1,
				const double mass2, const std::string cinfo,
				ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Uniform: public ISR_Channel_Base {
  protected:
    double m_exponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Simple_Pole_Uniform(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info,const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Central: public ISR_Channel_Base {
  protected:
    double m_exponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Simple_Pole_Central(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info,const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Forward: public ISR_Channel_Base {
  protected :
    double m_sexponent, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Simple_Pole_Forward(const double sexponent,const double yexponent,
			const std::string cinfo,ATOOLS::Integration_Info *info,
			const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Backward: public ISR_Channel_Base {
  protected :
    double m_sexponent, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Simple_Pole_Backward(const double sexponent,const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info,
			 const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };
}
#endif
