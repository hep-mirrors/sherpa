#ifndef LL_KPerp_H
#define LL_KPerp_H

#include "Single_Channel.H"

namespace PHASIC {

  class LL_KPerp: public Single_Channel {
  private:
    
    double m_beta;

    Info_Key m_kp1key, m_kp2key;

  public :

    // constructor
    LL_KPerp(double beta,const std::string cinfo,Integration_Info *info);

    // member functions
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode);
    void GeneratePoint(Info_Key &spkey,Info_Key &ykey,const int mode);
    void GenerateWeight(const int mode); 

    void CalculateLimits(Info_Key &spkey,Info_Key &ykey); 

  };// end of class LL_Kperp

}// end of namespace PHASIC

#endif
