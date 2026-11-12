#ifndef ATOOLS_Phys_Hadron_Collision_H
#define ATOOLS_Phys_Hadron_Collision_H

#include "ATOOLS/Math/Poincare.H"
#include<array>

namespace ATOOLS {
  class Particle;
  
  class Hadron_Collision {
  private:
    std::array<Particle *,2> m_particles;
    std::array<Vec4D,2>      m_momentaCMS, m_positionsCMS, m_velocitiesCMS;
    IntoOrientedCoM          m_Lorentz;
    Poincare_Sequence        m_InvLorentz;
    Vec4D                    m_momentumLab, m_positionLab, m_positionCMS;
    double                   m_s, m_tCMS, m_tLab, m_b2;
    bool                     m_on;
  public:
    Hadron_Collision(Particle * A,Particle * B);
    ~Hadron_Collision() {}

    Particle * InPart(const size_t pos)    const { return m_particles[pos]; }
    const Vec4D & MomCMS(const size_t pos) const { return m_momentaCMS[pos]; }

    const bool          On()    const { return m_on; }
    const double        T_Lab() const { return m_tLab; }
    const double        B2()    const { return m_b2; }
    const double        Ecms2() const { return m_s; }
    const ATOOLS::Vec4D X_Lab() const { return m_positionLab; }

    const Poincare_Sequence & InvBoost() const { return m_InvLorentz; }
  };
}

#endif
