#ifndef PHASIC_Channels_ISR_Channel_Base_H
#define PHASIC_Channels_ISR_Channel_Base_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/Info_Key.H"

namespace PHASIC {
  class ISR_Channel_Base: public Single_Channel {
  protected:
    Vegas  * p_vegas;
    double * p_rans;
    
    inline double SelectS(const double &s1,const double &s2) {
      return (s2>0.0?s2:s1);
    }
  public:
    ISR_Channel_Base(ATOOLS::Integration_Info *info) {}
    ~ISR_Channel_Base() {
      if (p_vegas) delete p_vegas;
      if (p_rans)  delete p_rans;
    }
    
    inline const std::string & ChID() const { return m_name; }
    inline void MPISync()                   { p_vegas->MPISync(); }
    inline void Optimize()                  { p_vegas->EndOptimize(); }
    inline void EndOptimize()               { p_vegas->EndOptimize(); }
    inline void WriteOut(std::string pId)   { p_vegas->WriteOut(pId); }
    inline void ReadIn(std::string pId)     { p_vegas->ReadIn(pId); }
    inline void AddPoint(double value) {
      Single_Channel::AddPoint(value);
      p_vegas->AddPoint(value,p_rans);
    }
  };
}
#endif
