#ifndef PHASIC_Channels_Leading_Log_Channels_H
#define PHASIC_Channels_Leading_Log_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"


namespace PHASIC {
  class Leading_Log_Uniform: public ISR_Channel_Base {
  protected:
    double m_beta, m_factor;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Leading_Log_Uniform(const double beta,const double factor,
			const std::string cinfo,ATOOLS::Integration_Info *info,
			const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };
  
  ////////////////////////////////////////////////////////////////////////

  class Leading_Log_Forward: public ISR_Channel_Base {
  protected :
    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Leading_Log_Forward(const double beta,const double factor,const double yexponent,
			const std::string cinfo,ATOOLS::Integration_Info *info,
			const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };
  
  ////////////////////////////////////////////////////////////////////////

  class Leading_Log_Backward: public ISR_Channel_Base {
  protected :
    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Leading_Log_Backward(const double beta,const double factor,const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info,
			 const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Leading_Log_Central: public ISR_Channel_Base {
  protected:
    double m_beta, m_factor;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Leading_Log_Central(const double beta,const double factor,
			const std::string cinfo,ATOOLS::Integration_Info *info,
			const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  ////////////////////////////////////////////////////////////////////////

  class Leading_Log_YFS: public ISR_Channel_Base {
  protected:
    double m_beta, m_factor;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_sgridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Leading_Log_YFS(const double beta,const double factor,
      const std::string cinfo,ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };
}
#endif
