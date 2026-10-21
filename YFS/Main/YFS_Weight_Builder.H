#ifndef YFS_Main_YFS_Weight_Builder_H
#define YFS_Main_YFS_Weight_Builder_H

#include "ATOOLS/Phys/Weights.H"

namespace YFS {

  class YFS_Handler;

  //! CalculateBeta through BuildNamedWeights and their ACRAIC/NLO probes,
  //! moved out of YFS_Handler. A friend of YFS_Handler, not a standalone
  //! calculation: it reads and writes YFS_Handler's own state directly.
  class YFS_Weight_Builder {
  public:
    explicit YFS_Weight_Builder(YFS_Handler &h) : m_h(h) {}

    void CalculateBeta();
    void GenerateWeight();

  private:
    YFS_Handler &m_h;

    void AcraicCompare();
    void InitNLO();
    double CalculateNLO();
    double AcraicRealVirtualRemainder();
    double AcraicRealVirtualHelicityRemainder();
    double AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const;
    void Beta02DumpProbe(double v) const;
    void BuildNamedWeights(double w_lo, double w_full);
    bool BuildLadderWeights(ATOOLS::Weights &w);

    void ProbeSign(double w_full);
    void ProbeLOFinal(double w_lo, double w_full, bool acraic_nom,
                      double corr_acraic, double corr_eex);
    void ProbeFinalWeights(double corr_eex, double corr_acraic, double wif,
                          bool acraico2, double corr_acraic_o1,
                          double corr_acraic_v0, double corr_acraic_v1);
    void ProbeLO();
    void ProbeLOPre(double pre_yfsw, double pre_ff, double pre_v, double pre_born);
    void ProbeNNLOEvent();
    void ProbeRealStab();
    void ProbeReference(bool haveacraic);
    void ProbeReferenceWeight(double w_lo);
    void ProbeAcraicCompare(bool haveacraic);
    void ProbeAcraicWeight();
    void ProbeGolden(double r0, double r1);
  };

}

#endif
