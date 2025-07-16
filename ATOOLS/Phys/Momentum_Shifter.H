#ifndef ATOOLS_Phys_Momentum_Shifter_H
#define ATOOLS_Phys_Momentum_Shifter_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"

#include <set>

namespace ATOOLS {

  struct ms {

    enum error_code {
      no_error        =   0,
      shift_error     =   1,
      direction_error =   2,
      sperp_error     =   4,
      momenta_error   =   8,
      boost_error_1   =  16,
      boost_error_2   =  32,
      scale_error     =  64,
      setup_error     = 128
    };

  };// end of struct ms

  std::ostream &operator<<(std::ostream &ostr,const ms::error_code code);

  class Momentum_Shifter {
  private:

    Vec4D m_shift, m_direction;
    Vec4D m_pold[3], m_pnew[3], m_pperp[3];

    Poincare  m_oldcms, m_newcms, m_rotate;
    Particle *p_initial[2];

    std::set<Particle*> m_boosted;

    size_t m_maxdepth;
    double m_sp[3];
    bool   m_setshift, m_setdirection, m_setsp[3];
    bool   m_initboost, m_initscale;

    bool CalculateShift();
    bool DetermineDirection();
    bool CalculateSPerp();
    bool ConstructMomenta();

    bool Boost(Particle *const particle,const size_t catcher);
    bool BoostBack(Particle *const particle,const size_t catcher);

  public:

    // constructor
    Momentum_Shifter(Particle *const initial1,Particle *const initial2);

    // member functions
    double Lambda2(double sp,double sp1,double sp2);

    bool Boost(Particle *const particle);
    bool BoostBack(Particle *const particle);
 
    ms::error_code Boost();
    ms::error_code BoostBack();

    ms::error_code Scale();

    void Reset();

    // inline functions
    inline void SetMaxDepth(const size_t depth)    { m_maxdepth=depth;    }

    inline void SetSPerp(const double sp,const size_t i) { m_sp[i]=sp; m_setsp[i]=true; }

    inline void SetShift(const Vec4D &shift)   { m_shift=shift; m_setshift=true;       }
    inline void SetDirection(const Vec4D &dir) { m_direction=dir; m_setdirection=true; }

    inline size_t MaxDepth() const { return m_maxdepth; }

    inline double SPerp(const size_t i) const { return m_sp[i]; }

    inline const Vec4D &Shift() const         { return m_shift;     }
    inline const Vec4D &Direction() const     { return m_direction; }

  };// end of class Momentum_Shifter

}// end of namespace ATOOLS

#endif
