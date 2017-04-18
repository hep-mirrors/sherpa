#ifndef AHADIC_Tools_Splitting_Tools_H
#define AHADIC_Tools_Splitting_Tools_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Vec.H"
#include "ATOOLS/Math/Poincare.H"

namespace AHADIC {
  class Splitting_Tools {
  private:
    static ATOOLS::Vec4D s_AxisP(1.,0.,0.,1.);
    static ATOOLS::Vec4D s_AxisM(1.,0.,0.,-1.);

    ATOOLS::Poincare m_boost, m_rotat;
    ATOOLS::Vec4D    m_splitmom,m_spectmom;
  public:
    Splitting_Tools();

    ConstructLightCone(ATOOLS::Particle * part1,ATOOLS::Particle * part2);
  };
  
  inline double Splitting_Tools::
  Lambda(const double & ma2,const double & mb2,const double & mc2) {
    return sqrt(ATOOLS::sqr(ma2-mb2-mc2)-4.*mb2*mc2)/(2.*ma2);
  };
}

#endif
