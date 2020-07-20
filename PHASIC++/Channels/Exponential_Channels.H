#ifndef PHASIC_Channels_Exponential_Channels_H
#define PHASIC_Channels_Exponential_Channels_H

#include "PHASIC++/Channels/ISR_Channel_Base.H"

namespace PHASIC {
  class Exponential_RelicDensity: public ISR_Channel_Base {
  protected:
    double m_exponent, m_mass[2];
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_sgridkey;
  public:
    Exponential_RelicDensity(const double exponent, const double mass1,const double mass2,
			     const std::string cinfo, ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };

  class Exponential_DM_Annihilation: public ISR_Channel_Base {
  protected:
    double m_exponent, m_mass[2];
    bool   m_zchannel;
    ATOOLS::Info_Key m_spkey, m_sgridkey, m_xkey, m_xgridkey, m_cosxikey, m_cosgridkey;
  public:
    Exponential_DM_Annihilation(const double exponent, const double mass1,const double mass2,
			     const std::string cinfo, ATOOLS::Integration_Info *info);
    void GeneratePoint(const double *rans);
    void GenerateWeight(const int & mode=0);
  };
} // end of namespace
#endif
