#ifndef PHASIC_Channels_Threshold_H
#define PHASIC_Channels_Threshold_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {

  class Threshold_Uniform: public Single_Channel {
  protected:

    double m_mass;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Threshold_Uniform(const double mass,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Threshold_Uniform

  class Threshold_Forward: public Single_Channel {
  protected :

    double m_mass, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Threshold_Forward(const double mass,const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Threshold_Forward

  class Threshold_Backward: public Single_Channel {
  protected :

    double m_mass, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Threshold_Backward(const double mass,const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Threshold_Backward

  class Threshold_Central: public Single_Channel {
  protected:

    double m_mass;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Threshold_Central(const double mass,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Threshold_Central

}// end of namespace PHASIC

#endif
