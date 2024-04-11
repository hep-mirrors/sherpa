#ifndef SHRIMPS_Ladders_Ladder_Generator_QT_H
#define SHRIMPS_Ladders_Ladder_Generator_QT_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_QT : public Ladder_Generator_Base {
  private:
    double        m_seff;
    double        m_sact, m_y[2];
    ATOOLS::Vec4D m_qini[2], m_cms, m_q[2];
    void CreateInitialMomenta();
    bool CreateLadder();
    void CreatePropagators();
    bool YStep();
    void MergeSinglets();
    void CalculateWeight();
  public:
    Ladder_Generator_QT();
    ~Ladder_Generator_QT() {}

    Ladder * operator()(const ATOOLS::Vec4D & pos);
  };
}
#endif
