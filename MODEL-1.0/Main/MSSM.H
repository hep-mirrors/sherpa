#ifndef MSSM_h
#define MSSM_h

#include "Model_Base.H"

namespace MODEL {
  class MSSM : public Model_Base {
  private :
    std::string m_generator;
    std::string m_benchmark,m_scenario;
    int         m_model;

    double      m_MW,m_MZ,m_MT;

    bool        m_scales,m_unification;
    int         m_sgnmu;
    double      m_tanbeta,m_A0,m_M0,m_M12;
    double      m_SoftMg[3],m_SoftM[3][5],m_SoftA[3],m_SoftH[2];
    double      m_SUSYScale;

    int         m_NMessengers;
    double      m_LambdaM,m_MScale,m_Gravitino;
    double      m_Gaugino,m_DeltaHu2,m_DeltaHd2;
    double      m_DeltaY,m_N5_1,m_N5_2,m_N5_3;

    double      m_M32;

    void InitializeBenchmarkPoint() {};
    void ReadInFile();
  public :
    MSSM(std::string,std::string);
    ~MSSM()                     {}
    bool   RunSpectrumGenerator();
    void   FillMasses()           {}
    void   FillWidths()           {}
  };
}

#endif


