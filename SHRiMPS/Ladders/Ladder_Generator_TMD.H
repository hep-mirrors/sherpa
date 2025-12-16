#ifndef SHRIMPS_Ladders_Ladder_Generator_TMD_H
#define SHRIMPS_Ladders_Ladder_Generator_TMD_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_TMD.H"
#include "REMNANTS/Tools/Pseudo_TMD.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_TMD : public Ladder_Generator_Base {
  private:
    double        m_S, m_Ecms, m_pt02;
    double        m_xprev[2], m_xtmdprev[2], m_thetaprev[2];
    double        m_ytest, m_qt2test, m_kttest, m_xtest, m_xtmdtest;
    size_t        m_ngluons;
    ATOOLS::Vec4D m_qtvectest, m_ktvectest, m_ktsum;
    
    REMNANTS::Pseudo_TMD * p_tmd[2];
    Sigma_TMD              m_sigmaTMD;
    
    void   CalculateWeight() {}
    void   InitialiseLadder(const ATOOLS::Vec4D & pos);
    void   FillZeroEmissionLadder();
    void   AddEmissions();
    bool   TrialEmission(const size_t beam);
    double EmissionWeight(const size_t beam);
    double AbsorptionWeight(const size_t beam);
    void   AddInitialState();
    void   AddPropagators();
    void   SelectPropagatorColours();
  public:
    Ladder_Generator_TMD();
    ~Ladder_Generator_TMD();

    void Initialise(PDF::ISR_Handler *const isr);
    Ladder * operator()(const ATOOLS::Vec4D & pos);
  };
}
#endif
