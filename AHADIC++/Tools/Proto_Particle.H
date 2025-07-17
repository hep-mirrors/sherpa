#ifndef AHADIC_Tools_Proto_Particle_H
#define AHADIC_Tools_Proto_Particle_H

#include "ATOOLS/Phys/Particle.H"
#include <iostream>

namespace AHADIC {
  class Proto_Particle {
  private:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_momentum, m_xprod, m_xdec;
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

    inline ATOOLS::Flavour Flavour()    const { return m_flav;     }
    inline ATOOLS::Vec4D   Momentum()   const { return m_momentum; }
    inline ATOOLS::Vec4D   XProd()      const { return m_xprod;    }
    inline ATOOLS::Vec4D   XDec()       const { return m_xdec;     }
    inline size_t          Generation() const { return m_gen; }
    inline double          KT2_Max()    const { return m_kt2max; }
    inline bool            IsLeading()  const { return m_isleading;  }
    inline bool            IsBeam()     const { return m_isbeam;  }

    inline void SetFlavour(const ATOOLS::Flavour & flav) { m_flav      = flav; }
    inline void SetMomentum(const ATOOLS::Vec4D & mom)   { m_momentum  = mom; }
    inline void SetXProd(const ATOOLS::Vec4D & xprod)    { m_xprod     = xprod; }
    inline void SetXDec(const ATOOLS::Vec4D & xdec)      { m_xdec      = xdec; }
    inline void SetGeneration(const size_t & gen)        { m_gen       = gen; }
    inline void SetKT2_Max(const double & kt2max)        { m_kt2max    = kt2max; }
    inline void SetLeading(const bool & lead)            { m_isleading = lead; }
    inline void SetBeam(const bool & beam)               { m_isbeam    = beam; }


    static std::set<Proto_Particle *> s_protos;
    static void Reset();

    ATOOLS::Particle * operator()();
    const ATOOLS::Vec4D Velocity() const;
  };

  std::ostream & operator<<(std::ostream &,const Proto_Particle &);
}
#endif
