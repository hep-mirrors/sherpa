#ifndef BEAM_Main_Kinematics_Base_H
#define BEAM_Main_Kinematics_Base_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include <map>
#include <array>

namespace BEAM {
  class Kinematics_Base {
  protected:
    std::array<Beam_Base*, 2> p_beams;
    bool        m_on;
    std::string m_keyid;
    std::array<double, 2> m_m, m_m2, m_exponent;
    /*!
      The energy fractions the outgoing bunches have w.r.t. the corresponding
      incoming beams.
    */
    std::array<double, 2> m_x;
    double m_smin, m_smax, m_sminPS, m_smaxPS, m_S;

    ATOOLS::Vec4D    m_Plab;
    /*!
      A boost from the c.m. system of the incoming beams to the c.m. system of the
      outgoing bunches.
    */
    ATOOLS::Poincare m_CMSBoost;
    /*!
      A boost from the lab system of the incoming beams to their c.m. system.
    */
    ATOOLS::Poincare m_LabBoost;
    virtual void InitIntegration()=0;
  public:
    Kinematics_Base(std::array<Beam_Base *, 2> beams);
    virtual ~Kinematics_Base();

    virtual bool operator()(ATOOLS::Vec4D_Vector& moms)           = 0;
    virtual void AssignKeys(ATOOLS::Integration_Info *const info) = 0;
    virtual void SetLimits()                                      = 0;
    virtual inline void SetSprimeMin(const double & smin) { m_sminPS = ATOOLS::Max(smin, m_sminPS); }
    virtual inline void SetSprimeMax(const double & smax) { m_smaxPS = smax; }
    virtual inline const double & SprimeMin() const { return ATOOLS::Max(m_smin, m_sminPS); }
    virtual inline const double & SprimeMax() const { return ATOOLS::Min(m_smax, m_smaxPS); }
    virtual inline double         Sprime()          { return m_S; }
    virtual inline double         Y()               { return 0.; }

    virtual inline std::array<Beam_Base*, 2> GetBeams() const { return p_beams; }

    virtual inline bool On()                           const { return m_on; }
    virtual inline const std::string & KeyId()         const { return m_keyid; }
    /**
     * @brief Gets the mass of one of the incoming particles.
     * @param i Which particle to get
     * @return the mass
     */
    virtual inline const double & m(const size_t & i)  const { return m_m[i]; }
    /**
     * @brief Gets the square of the mass of one of the incoming particles.
     * @param i Which particle to get
     * @return the mass
     */
    virtual inline const double & m2(const size_t & i) const { return m_m2[i]; }
    /*!
      Characteristic exponents used for the integration. This is mainly a service for better
      performance of the phase space integration.
    */
    virtual inline const double & Exponent(const size_t & i) const {
      return m_exponent[i];
    }
    virtual inline void BoostInCMS(ATOOLS::Vec4D & p) { m_CMSBoost.Boost(p); }
    virtual inline void BoostInCMS(ATOOLS::Vec4D_Vector& p) {
      for (auto& mom : p) m_CMSBoost.Boost(mom);
    }
    virtual inline void BoostInLab(ATOOLS::Vec4D & p) { m_LabBoost.Boost(p); }
    virtual inline void BoostInLab(ATOOLS::Vec4D_Vector& p) {
      for (auto& mom : p) m_LabBoost.Boost(mom);
    }
    virtual inline void BoostBackCMS(ATOOLS::Vec4D & p) { m_CMSBoost.BoostBack(p); }
    virtual inline void BoostBackCMS(ATOOLS::Vec4D_Vector& p) {
      for (auto& mom : p) m_CMSBoost.BoostBack(mom);
    }
    virtual inline void BoostBackLab(ATOOLS::Vec4D & p) { m_LabBoost.BoostBack(p); }
    virtual inline void BoostBackLab(ATOOLS::Vec4D_Vector& p) {
      for (auto& mom : p) m_LabBoost.BoostBack(mom);
    }
  };
}
#endif
