#ifndef SHRIMPS_EVENT_GENERATOR_BASE_H
#define SHRIMPS_EVENT_GENERATOR_BASE_H

#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace SHRIMPS {
  class Event_Generator_Base {
  protected:
    Omega_ik * p_eikonal;
    double     m_xsec, m_smin;
  public:
    Event_Generator_Base(Sigma_Base * sigma):
      p_eikonal(NULL), m_xsec(sigma->XSec()), m_smin(0.) {} 
    ~Event_Generator_Base() {};
    virtual void       Reset() {}
    virtual int        GenerateEvent(ATOOLS::Blob_List * blobs,const bool & flag=false) = 0;
    virtual Omega_ik * GetEikonal()      const { return p_eikonal;}
    virtual double     XSec()            const { return m_xsec; }
    virtual double     Smin()            const { return m_smin;}
    virtual bool       IsLastRescatter() const { return false; }
    virtual double     TMax()            const { return 0.; }
    virtual int        NLadders()        const { return 1; }
  };
}
#endif
