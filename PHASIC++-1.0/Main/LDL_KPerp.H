#ifndef LDL_KPerp_H
#define LDL_KPerp_H

#include "Single_Channel.H"

namespace PHASIC {

  class LDL_KPerp: public Single_Channel {
  private:
    
    double m_lambda2, m_rancut[2], m_integral[2], m_firstintegral[2];
    bool   m_calculated;

    ATOOLS::Info_Key m_kp1key, m_kp2key;

  public :

    // constructor
    LDL_KPerp(const double lambda2,const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode);
    void GenerateWeight(const int mode); 

    void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey); 
    void CalculateLimits(); 

  };// end of class LDL_Kperp

}// end of namespace PHASIC

#endif
