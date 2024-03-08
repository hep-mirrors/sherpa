#ifndef SHRIMPS_Ladders_Ladder_Generator_Eik_H
#define SHRIMPS_Ladders_Ladder_Generator_Eik_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_Eik : public Ladder_Generator_Base {
  private:
    ATOOLS::Vec4D   m_lastk;
    void   FillGluons();
    void   SelectPropagatorColours();
    bool   SelectPropagatorQTs();
    bool   SelectPropagatorQT(const size_t dir,T_Prop & prop);
    void   CalculateWeight();
  public:
    Ladder_Generator_Eik();
    ~Ladder_Generator_Eik();
    
    Ladder * operator()(const ATOOLS::Vec4D & pos);
  };
}
#endif
