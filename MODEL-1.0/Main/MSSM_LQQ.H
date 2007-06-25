#ifndef MSSM_LQQ_h
#define MSSM_LQQ_h

#include "Model_Base.H"

namespace MODEL {
  class MSSM_LQQ : public Model_Base {
  private :
    MSSM * p_MSSM;

    std::string m_generator;
    std::string m_benchmark,m_scenario;

    bool        m_scales,m_unification;

    void InitializeBenchmarkPoint() {};
    void ReadInFile();
  public :
    MSSM_LQQ(std::string,std::string);
    ~MSSM_LQQ()                     {}
    bool RunSpectrumGenerator();
    void FillMasses()               {}
    void FillWidths()               {}
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
}

#endif


