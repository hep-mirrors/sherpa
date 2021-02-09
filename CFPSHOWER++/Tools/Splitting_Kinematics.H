#ifndef CFPSHOWER_Tools_Splitting_Kinematics_H
#define CFPSHOWER_Tools_Splitting_Kinematics_H

#include "ATOOLS/Math/Poincare.H"

namespace CFPSHOWER {
  class Splitting;
  
  class Splitting_Kinematics_Base {
  public:
    Splitting_Kinematics_Base(Splitting & split);
    ATOOLS::Vec4D m_psplit, m_pspect, m_pboth;
    double        m_Q2, m_Q, m_m2[4], m_m[4];
  }
}

#endif

