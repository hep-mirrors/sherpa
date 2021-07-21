#ifndef BEAM_Main_Kinematics_Base_H
#define BEAM_Main_Kinematics_Base_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include <map>

namespace BEAM {
  class Kinematics_Base {
  protected:
    Beam_Base * p_beams[2];
    bool        m_on;
    std::string m_keyid;
    double      m_m[2], m_m2[2], m_x[2], m_S;
    double      m_smin, m_smax, m_sminPS, m_smaxPS, m_exponent[2];

    ATOOLS::Vec4D    m_Plab;
    ATOOLS::Poincare m_CMSBoost, m_LABBoost;
    virtual void InitIntegration()=0;
  public:
    Kinematics_Base(Beam_Base * p_beams[2]);
    ~Kinematics_Base();
    
    virtual bool operator()(ATOOLS::Vec4D * moms)                 = 0;
    virtual void AssignKeys(ATOOLS::Integration_Info *const info) = 0;
    virtual void SetLimits()                                      = 0;
    virtual inline void SetSprimeMin(const double & smin) { m_sminPS = smin; }
    virtual inline void SetSprimeMax(const double & smax) { m_smaxPS = smax; }
    virtual inline const double & SprimeMin() const { return ATOOLS::Max(m_smin, m_sminPS); }
    virtual inline const double & SprimeMax() const { return ATOOLS::Min(m_smax, m_smaxPS); }
    virtual inline const double   Sprime()          { return m_S; }
    virtual inline const double   Y()               { return 0.; }
    
    virtual inline Beam_Base ** GetBeams() { return p_beams; }

    virtual inline const bool     On()                 const { return m_on; }
    virtual inline const std::string & KeyId()         const { return m_keyid; }
    virtual inline const double & m(const size_t & i)  const { return m_m[i]; }
    virtual inline const double & m2(const size_t & i) const { return m_m2[i]; }
    virtual inline const double & Exponent(const size_t & i) const {
      return m_exponent[i];
    }

    virtual inline void BoostInCMS(ATOOLS::Vec4D* p,const size_t & n) {
      for (size_t i=0; i<n; ++i) m_CMSBoost.Boost(p[i]);
    }
    
    virtual inline void BoostInLab(ATOOLS::Vec4D* p,const size_t & n) {
      for (size_t i=0; i<n; ++i) m_CMSBoost.BoostBack(p[i]);
    }
  };
}
#endif
