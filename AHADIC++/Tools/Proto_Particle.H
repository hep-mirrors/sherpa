#ifndef AHADIC_Tools_Proto_Particle_H
#define AHADIC_Tools_Proto_Particle_H

#include "ATOOLS/Phys/Particle.H"
#include <iostream>

namespace AHADIC {
  class Proto_Particle {
  private:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_momentum;
    size_t          m_gen;
    double          m_kt2max;
    bool            m_isleading, m_isbeam;
  public:
    Proto_Particle(const Proto_Particle & proto);
    Proto_Particle(const ATOOLS::Particle & part);
    Proto_Particle(const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_none),
		   const ATOOLS::Vec4D & mom=ATOOLS::Vec4D(0.,0.,0.,0.),
		   bool leading=false,bool beam=false);
    ~Proto_Particle();

    ATOOLS::Flavour Flavour()    const { return m_flav;     }
    ATOOLS::Vec4D   Momentum()   const { return m_momentum; }
    size_t          Generation() const { return m_gen; }
    double          KT2_Max()    const { return m_kt2max; }
    bool            IsLeading()  const { return m_isleading;  }
    bool            IsBeam()     const { return m_isbeam;  }

    void SetFlavour(const ATOOLS::Flavour & flav) { m_flav      = flav; }
    void SetMomentum(const ATOOLS::Vec4D & mom)   { m_momentum  = mom; }
    void SetGeneration(const size_t & gen)        { m_gen       = gen; }
    void SetKT2_Max(const double & kt2max)        { m_kt2max    = kt2max; }
    void SetLeading(const bool & lead)            { m_isleading = lead; }
    void SetBeam(const bool & beam)               { m_isbeam    = beam; }


    static std::set<Proto_Particle *> s_protos;
    static void Reset();

    ATOOLS::Particle * operator()();
  };

  std::ostream & operator<<(std::ostream &,const Proto_Particle &);
}
#endif
