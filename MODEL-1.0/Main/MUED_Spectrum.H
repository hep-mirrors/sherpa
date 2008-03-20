#ifndef MUED_Spectrum_h
#define MUED_Spectrum_h

#include "Spectrum_Generator_Base.H"
#include "Model_Base.H"
#include "Data_Reader.H"
#include "Decay_Table.H"

namespace MODEL {
  class MUED_Spectrum : public Spectrum_Generator_Base {
  private:
    int    m_generations;
    double m_invR, m_Lambda, m_Lambda2, m_Riemann3;
    double m_sin2thetaW, m_cos2thetaW, m_v;

    void Neutral_KK_Bosons();
    void Charged_KK_Bosons();
    void Neutral_KK_Scalars();
    void Charged_KK_Scalars();
    void LR_KK_Fermions();
  public:
    MUED_Spectrum(ATOOLS::Data_Reader *,Model_Base *);
    ~MUED_Spectrum();
    void Run(std::string);
    void FillMasses();
    void FillDecays()                     { }
    bool FillDecay(ATOOLS::Decay_Table *) { return true; }
  };
}

#endif


