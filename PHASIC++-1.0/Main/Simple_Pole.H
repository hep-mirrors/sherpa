#ifndef Simple_Pole_H
#define Simple_Pole_H

#include "Single_Channel.H"

namespace PHASIC {

  class Simple_Pole_Uniform: public Single_Channel {
  protected:

    double m_exponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Simple_Pole_Uniform(const double exponent,const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Simple_Pole_Uniform

  class Simple_Pole_Forward: public Single_Channel {
  protected :

    double m_sexponent, m_yexponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Simple_Pole_Forward(const double sexponent,const double yexponent,
		      const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Simple_Pole_Forward

  class Simple_Pole_Backward: public Single_Channel {
  protected :

    double m_sexponent, m_yexponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Simple_Pole_Backward(const double sexponent,const double yexponent,
		       const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Simple_Pole_Backward

  class Simple_Pole_Central: public Single_Channel {
  protected:

    double m_exponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Simple_Pole_Central(const double exponent,const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Simple_Pole_Central

}// end of namespace PHASIC

#endif
