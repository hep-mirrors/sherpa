#ifndef Leading_Log_Z_H
#define Leading_Log_Z_H

#include "Single_Channel.H"
#include "Vegas.H"

namespace PHASIC {

  class Z_Channel_Base: public Single_Channel {
  protected:

    double m_beta, m_cut; 

    Vegas *p_vegas;

    double f(ATOOLS::Info_Key &key,double t);
    double I(ATOOLS::Info_Key &key,double t);
    double t(ATOOLS::Info_Key &key,double I);

  public:

    // constructor
    Z_Channel_Base();
    
    // destructor
    ~Z_Channel_Base();

    // member functions
    std::string ChID();

    void Optimize();
    void EndOptimize();
    void WriteOut(std::string pId);
    void ReadIn(std::string pId);

  };// end of class Z_Channel_Base

  class Leading_Log_Z_QQ: public Z_Channel_Base {
  protected:

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_QQ(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double value);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_QQ

  class Leading_Log_Z_QG: public Z_Channel_Base {
  protected:

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_QG(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double value);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_QG

  class Leading_Log_Z_GQ: public Z_Channel_Base {
  protected:

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_GQ(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double value);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_GQ

  class Leading_Log_Z_GG: public Z_Channel_Base {
  protected:

    ATOOLS::Info_Key m_zkey[2];

  public:

    // constructor
    Leading_Log_Z_GG(const double beta,const double cut,
		     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double value);

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);

  };// end of class Leading_Log_Z_GG

}// end of namespace PHASIC

#endif
