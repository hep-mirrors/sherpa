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
    ATOOLS::Vec4D     m_pbeam[2], m_ktsum;
    double            m_kt2max, m_phi, m_yseed[4];
    size_t            m_emissions[3];
    colour_type::code m_cols[3];
    

    void   SeedLadder(const ATOOLS::Vec4D & pos);
    void   AddBeamPartons();
    void   FillIntervals();
    void   CompensateKTs();
    void   ConstructFSMomenta();
    void   FillPropagators();
    ATOOLS::Vec4D SelectKT(const double & y);
    void   CalculateWeight();
  public:
    Ladder_Generator_Seeded();
    ~Ladder_Generator_Seeded();
    
    Ladder * operator()(const ATOOLS::Vec4D & pos);    
  };
}
#endif
