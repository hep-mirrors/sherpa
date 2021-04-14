#ifndef SHRIMPS_Ladders_Ladder_Generator_KT_H
#define SHRIMPS_Ladders_Ladder_Generator_KT_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_KT : public Ladder_Generator_Base {
  private:
    Continued_PDF * p_pdf[2];
    double          m_Ebeam[2];
    double          m_S, m_maxsigma, m_ktmin, m_shatmin, m_shat;
    double          m_seff, m_kt2, m_weight;
    
    ATOOLS::Vec4D       m_q[2], m_k;
    TPropList::iterator m_pit[2];
    
    void   ScanPDFs();
    bool   FixInitialState();
    bool   FillForwardPartons();
    bool   FillEmissions();
    bool   MakeTrialLadder();
    bool   TrialEmission(double & y,const double & arg,const size_t dir);
    double AbsorptionWeight(const ATOOLS::Vec4D & k,const double & deltay);
    bool   Ordering(const size_t dir);
    void   AddEmission(const size_t dir,const double & y);
    bool   AdjustLastEmission(const bool dlast, const double & ylast,const double & y0);
    bool   UpdateInitialState();    
    bool   FinishKinematics();    
    bool   FixSimpleKinematics();
    void   SelectPropagatorColours();
    void   CalculateWeight();
    double KTMaxWeight();
  public:
    Ladder_Generator_KT();
    ~Ladder_Generator_KT() {}
    
    Ladder * operator()(const ATOOLS::Vec4D & pos);
    inline void SetRemnantHandler(Remnant_Handler * remnants) {
      for (size_t beam=0;beam<2;beam++) p_pdf[beam] = remnants->GetPDF(beam);
    }

  };
}
#endif
