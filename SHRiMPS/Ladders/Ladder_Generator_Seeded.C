//emission k_t's can become small, maybe such emissions need to be kicked out
#include "SHRiMPS/Ladders/Ladder_Generator_Seeded.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_Seeded::Ladder_Generator_Seeded() : Ladder_Generator_Base() {}

Ladder_Generator_Seeded::~Ladder_Generator_Seeded() {}

Ladder * Ladder_Generator_Seeded::operator()(const Vec4D & pos) {
  SeedLadder();
  return p_ladder;
}

void Ladder_Generator_Seeded::SeedLadder() {
}


double Ladder_Generator_Seeded::PDFFactor() {
  for (size_t beam=0;beam<2;beam++) {
    m_pdf[beam][0] = 0.;
    for (long int fl=-5;fl<6;fl++) {
      Flavour flav    = (fl==0?Flavour(kf_gluon):Flavour(fl));
      m_pdf[beam][0] += m_pdf[beam][6+fl] = m_partonic.PDF(beam,m_xseed[beam],m_pt2,flav);
    }
  }
  return m_pdf[0][0]*m_pdf[1][0];
}

void Ladder_Generator_Seeded::CalculateWeight() {
  m_weight  = 1.;
}

