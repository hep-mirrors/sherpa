#ifndef SHRIMPS_Ladders_Ladder_Generator_Seeded_H
#define SHRIMPS_Ladders_Ladder_Generator_Seeded_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_Seeded : public Ladder_Generator_Base {
  private:
    double m_pt2, m_yseed[2], m_xseed[2], m_pdf[2][12];

    void   SeedLadder();
    double PDFFactor();
    void   CalculateWeight();
  public:
    Ladder_Generator_Seeded();
    ~Ladder_Generator_Seeded();
    
    Ladder * operator()(const ATOOLS::Vec4D & pos);    
  };
}
#endif
