#ifndef PHASIC_Channels_LBS_Compton_Peak_H
#define PHASIC_Channels_LBS_Compton_Peak_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {

  class LBS_Compton_Peak_Uniform: public Single_Channel {
  protected:

    double m_exponent, m_pole;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    LBS_Compton_Peak_Uniform(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info);
    
    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class LBS_Compton_Peak_Uniform

  class LBS_Compton_Peak_Forward: public Single_Channel {
  protected :

    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    LBS_Compton_Peak_Forward(const double exponent,const double pole,const double yexponent,
			     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class LBS_Compton_Peak_Forward

  class LBS_Compton_Peak_Backward: public Single_Channel {
  protected :

    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    LBS_Compton_Peak_Backward(const double exponent,const double pole,const double yexponent,
			      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class LBS_Compton_Peak_Backward

  class LBS_Compton_Peak_Central: public Single_Channel {
  protected:

    double m_exponent, m_pole;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    LBS_Compton_Peak_Central(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class LBS_Compton_Peak_Central

}// end of namespace PHASIC

#endif
