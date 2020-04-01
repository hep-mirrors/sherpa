#ifndef PHASIC_Channels_Threshold_Channels_H
#define PHASIC_Channels_Threshold_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"
#include "PHASIC++/Channels/Vegas.H"


namespace PHASIC {
  class Threshold_Uniform: public ISR_Channel_Base {
  protected:
    double m_mass,m_sexp;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Threshold_Uniform(const double mass,const double sexp,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
    inline void   SetMass(const double mass) { m_mass=mass; }
    inline double Mass() const { return m_mass; }
  };

  ////////////////////////////////////////////////////////////////////////

  class Threshold_Forward: public ISR_Channel_Base {
  protected :
    double m_mass,m_sexp, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Threshold_Forward(const double mass,const double sexp,
		      const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
    inline void   SetMass(const double mass) { m_mass=mass; }
    inline const double Mass()         const { return m_mass; }
  };

  ////////////////////////////////////////////////////////////////////////

  class Threshold_Backward: public ISR_Channel_Base {
  protected :
    double m_mass,m_sexp, m_yexponent;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public :
    Threshold_Backward(const double mass,const double sexp,
		       const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info,
		       const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
    inline void   SetMass(const double mass) { m_mass=mass; }
    inline const double Mass()         const { return m_mass; }
  };

  ////////////////////////////////////////////////////////////////////////

  class Threshold_Central: public ISR_Channel_Base {
  protected:
    double m_mass,m_sexp;
    bool   m_zchannel;
    size_t m_mode;
    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    ATOOLS::Info_Key m_kp1key, m_kp2key;
  public:
    Threshold_Central(const double mass,const double sexp,
		      const std::string cinfo,ATOOLS::Integration_Info *info,
		      const size_t mode=0);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
    inline void   SetMass(const double mass) { m_mass=mass; }
    inline const double Mass()         const { return m_mass; }
  };

}
#endif
