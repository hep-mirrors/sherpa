#ifndef LL_KPerp_H
#define LL_KPerp_H

#include "Single_Channel.H"
#include "Vegas.H"

namespace PHASIC {

  class LL_KPerp: public Single_Channel {
  private:
    
    double m_beta, m_rancut[2], m_integral[2], m_firstintegral[2];
    bool   m_calculated;

    ATOOLS::Info_Key m_kpkey[2];

    Vegas *p_vegas;

    double f(ATOOLS::Info_Key &key,double t);
    double I(ATOOLS::Info_Key &key,double t);
    double t(ATOOLS::Info_Key &key,double I);

  public :

    // constructor
    LL_KPerp(double beta,const std::string cinfo,ATOOLS::Integration_Info *info);

    // destructor
    ~LL_KPerp();

    std::string ChID();

    void Optimize();
    void EndOptimize();
    void WriteOut(std::string pId);
    void ReadIn(std::string pId);

    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode);
    void GenerateWeight(const int mode); 

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey); 
    void AddPoint(double value);

  };// end of class LL_Kperp

}// end of namespace PHASIC

#endif
