#ifndef Leading_Log_H
#define Leading_Log_H

#include "Single_Channel.H"

namespace PHASIC {

  class Leading_Log_Uniform: public Single_Channel {
  protected:

    double m_beta, m_factor;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Leading_Log_Uniform(const double beta,const double factor,
			const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    bool CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Leading_Log_Uniform

  class Leading_Log_Forward: public Single_Channel {
  protected :

    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Leading_Log_Forward(const double beta,const double factor,const double yexponent,
			const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    bool CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Leading_Log_Forward

  class Leading_Log_Backward: public Single_Channel {
  protected :

    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public :

    // constructor
    Leading_Log_Backward(const double beta,const double factor,const double yexponent,
			 const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    bool CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Leading_Log_Backward

  class Leading_Log_Central: public Single_Channel {
  protected:

    double m_beta, m_factor;
    bool   m_zchannel;

    Info_Key m_spkey, m_xkey, m_ykey;

  public:

    // constructor
    Leading_Log_Central(const double beta,const double factor,
			const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    bool CalculateLimits(Info_Key &spkey,Info_Key &ykey);

  };// end of class Leading_Log_Central

}// end of namespace PHASIC

#endif
