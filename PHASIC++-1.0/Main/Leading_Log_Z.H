#ifndef Leading_Log_Z_H
#define Leading_Log_Z_H

#include "Single_Channel.H"

namespace PHASIC {

  class Leading_Log_Z_QQ: public Single_Channel {
  protected:

    double m_beta, m_cut;

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_QQ(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_QQ

  class Leading_Log_Z_QG: public Single_Channel {
  protected:

    double m_beta, m_cut;

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_QG(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_QG

  class Leading_Log_Z_GQ: public Single_Channel {
  protected:

    double m_beta, m_cut;

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_GQ(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_GQ

  class Leading_Log_Z_GG: public Single_Channel {
  protected:

    double m_beta, m_cut;

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_GG(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_GG

}// end of namespace PHASIC

#endif
