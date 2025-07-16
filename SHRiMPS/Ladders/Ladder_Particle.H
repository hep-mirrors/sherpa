#ifndef SHRIMPS_Ladders_Ladder_Particle_H
#define SHRIMPS_Ladders_Ladder_Particle_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vector.H"
#include <map>

namespace SHRIMPS {
  class Ladder_Particle {
  private:
    ATOOLS::Flavour    m_flav;
    ATOOLS::Vec4D      m_mom, m_pos;
    ATOOLS::Flow       m_flow;
    bool               m_marked, m_IS;
    int                m_beam;
  public:
    Ladder_Particle(const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_none),
		    const ATOOLS::Vec4D & mom=ATOOLS::Vec4D(0.,0.,0.,0.),
		    const ATOOLS::Vec4D & pos=ATOOLS::Vec4D(0.,0.,0.,0.));
    Ladder_Particle(const Ladder_Particle & part);
    ~Ladder_Particle() { }

    void SetFlow(const unsigned int & pos,const int & code=-1);
    ATOOLS::Particle * GetParticle();
    
    inline const ATOOLS::Flavour & Flavour()       const { return m_flav; }
    inline const ATOOLS::Vec4D & Momentum()        const { return m_mom; }
    inline const ATOOLS::Vec4D & Position()        const { return m_pos; }
    inline const int GetFlow(const size_t & pos)   const { return m_flow.Code(pos); }
    inline const bool & IsIS()                     const { return m_IS; }
    inline const bool & IsMarked()                 const { return m_marked; }
    inline const int & Beam()                      const { return m_beam; }
    
    inline void SetFlavour(const ATOOLS::Flavour & flav)    { m_flav   = flav; }
    inline void SetMomentum(const ATOOLS::Vec4D & mom)      { m_mom    = mom; }
    inline void SetPosition(const ATOOLS::Vec4D & pos)      { m_pos    = pos; }
    inline void SetIS(const bool & IS)                      { m_IS     = IS; }
    inline void SetMark(const bool & mark)                  { m_marked = mark; }
    inline void SetBeam(const int & beam)                   { m_beam   = beam; }
    inline void BoostBack(const ATOOLS::Poincare & booster) { booster.BoostBack(m_mom); }
    inline void Boost(const ATOOLS::Poincare & booster)     { booster.Boost(m_mom); }
  };
  typedef std::map<double,Ladder_Particle,std::greater<double> > LadderMap;
  std::ostream & operator<<(std::ostream & s, const Ladder_Particle & part);
  std::ostream & operator<<(std::ostream & s, const LadderMap & lmap);
}

#endif
