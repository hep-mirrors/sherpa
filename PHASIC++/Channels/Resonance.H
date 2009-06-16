#ifndef PHASIC_Channels_Resonance_H
#define PHASIC_Channels_Resonance_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {

  class Resonance_Uniform: public Single_Channel {
  protected:

    double m_mass, m_width;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Resonance_Uniform(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Resonance_Uniform

  class Resonance_Forward: public Single_Channel {
  protected :

    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Resonance_Forward(const double mass,const double width,const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Resonance_Forward

  class Resonance_Backward: public Single_Channel {
  protected :

    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Resonance_Backward(const double mass,const double width,const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Resonance_Backward

  class Resonance_Central: public Single_Channel {
  protected:

    double m_mass, m_width;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Resonance_Central(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

  };// end of class Resonance_Central

}// end of namespace PHASIC

#endif
