#ifndef AHADIC_Tools_Multiplet_Constructor_H
#define AHADIC_Tools_Multiplet_Constructor_H

#include "AHADIC++/Tools/Hadron_Multiplet.H"


namespace AHADIC {
  class Hadronisation_Reweighting;

  struct HadInfo {
    ATOOLS::Flavour flav;
    int         iso,exr,exl,fl3,fl2,fl1,spin2;
    bool        barrable;
    std::string multiname;
    double      spinwt;
    // per-variation multiplet and extra weights for the reweighting
    std::vector<double> multiwt, extrawt;

    friend std::ostream & operator<<(std::ostream & s,const HadInfo & info);
  };


  class Multiplet_Constructor {
  private:
    bool    m_test;
    Hadronisation_Reweighting * p_reweighting;
    size_t  m_n_soft_variations {0};
    std::vector<double> m_singletsuppression, m_etam, m_etapm;
    std::vector<double> m_cse, m_bse, m_bce, m_hcbe, m_hbbe, m_sbm;
    std::vector<double> m_mixing0, m_mixing1, m_mixing2;
    HadInfo m_info;
    
    std::map<std::string,std::set<ATOOLS::Flavour> > m_multiplets;
    Wave_Functions m_wavefunctions;
    
    void        CreateMultiplets();
    bool        ExtractInfo(const kf_code & kfc);
    std::string MultipletName();
    void        FillMultipletWeights();
    void        ScaleExtraWeights(const std::vector<double> & factors);
    bool        ComponentActive(const std::vector<double> & amps);
    std::vector<double> MixingAngles();

    bool ConstructWaveFunction();
    bool ConstructMesonWaveFunction();
    Wave_Function * TrivialMesonWaveFunction();
    Wave_Function * Pi0WaveFunction();
    Wave_Function * OctetMesonWaveFunction();
    Wave_Function * SingletMesonWaveFunction();
    bool ConstructBaryonWaveFunction();
    bool ConstructAntiBaryonWaveFunction(ATOOLS::Flavour & flav);
    Wave_Function * NeutronWaveFunction();
    Wave_Function * ProtonWaveFunction();
    Wave_Function * SigmaWaveFunction();
    Wave_Function * LambdaWaveFunction();
    Wave_Function * Lambda1WaveFunction();
    Wave_Function * SigmaHWaveFunction();
    Wave_Function * LambdaHWaveFunction();
    Wave_Function * DeltaPPWaveFunction();
    Wave_Function * DeltaPWaveFunction();
    Wave_Function * Delta0WaveFunction();
    Wave_Function * DecupletSigmaWaveFunction();
    Wave_Function * DecupletLambdaWaveFunction();
    Wave_Function * HeavyDecupletWaveFunction();
  public:
    Multiplet_Constructor(bool test, Hadronisation_Reweighting * reweighting);
    ~Multiplet_Constructor();

    Wave_Functions * GetWaveFunctions() { return &m_wavefunctions; }  
    void PrintWaveFunctions(bool checkonly=false);
    void PrintMultiplets();
  };
}

#endif
