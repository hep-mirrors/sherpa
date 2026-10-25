#ifndef YFS_Main_YFS_Weight_Builder_H
#define YFS_Main_YFS_Weight_Builder_H

#include "ATOOLS/Phys/Weights.H"
#include <vector>
#include <string>

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

    //! The ACRAIC named-column pieces GenerateWeight builds, threaded
    //! through ApplyExternalVirtual/ApplyHelicityVirtual and friends.
    struct AcraicColumns {
      double corr_acraic = 0.;
      double corr_acraic_o1 = 0., corr_acraic_o1v1 = 0., corr_acraic_o1s = 0.;
      double corr_acraic_v0 = 0., corr_acraic_v1 = 0., corr_acraic_v2 = 0.;
      double corr_acraic_o0 = 0., corr_acraic_o0_ifcut = 0.;
      double corr_acraic_fo = 0.;
      double corr_acraic_nob02 = 0., corr_acraic_novp = 0.;
      double tnp_vhalf = 0., vext_fo = 0., dyif_fo = 0., dv2_acraic = 0.;
      double b02_acraic = 0., b02_v0 = 0., b02_v1 = 0.;
      double vp_acraic = 0., vp_v0 = 0., vp_v1 = 0.;
      std::vector<double> tnp1, tnp2;
      bool acraico2 = false, acraicb02 = false, acraicvp = false;
    };
    //! ACRAIC_Virtual: external/automatic - the squared-level provider's v on A_0
    void ApplyExternalVirtual(AcraicColumns &c);
    //! ACRAIC_Virtual: helicity - the same insertion, per helicity
    void ApplyHelicityVirtual(AcraicColumns &c);
    //! beta_0^(2) and Dyson VP, added after either virtual above; TNP templates
    void FinalizeAcraicBeta02VP(AcraicColumns &c);
    //! ACRAIC: IF_FORMFACTOR_CUT
    void ApplyIFFormFactorCut(AcraicColumns &c);
    //! ACRAIC: FIXED_ORDER_STREAM
    void ComputeFixedOrderColumn(AcraicColumns &c);

    void AcraicCompare();
    void InitNLO();
    double CalculateNLO();
    double AcraicRealVirtualRemainder();
    double AcraicRealVirtualHelicityRemainder();
    double AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const;
    void BuildNamedWeights(double w_lo, double w_full);
    //! the YFS.* column names this configuration emits, and wyfs/m_wnames defaults
    void RegisterColumnNames(std::vector<std::string> &names,
                             std::vector<std::string> &tnpnames,
                             ATOOLS::Weights &wyfs);
    void FillAcraicColumns(ATOOLS::Weights &wyfs);
    void FillTnpColumns(ATOOLS::Weights &wyfs, const std::vector<std::string> &tnpnames);
    //! YFS: DIFF_BUDGET_STREAMS
    void FillDiffBudgetColumns(ATOOLS::Weights &wyfs);
    //! NoCoulombHO, ACRAIC_NO_B02/VP, NLO_EW_NO_B02/VP
    void FillWithoutColumns(ATOOLS::Weights &wyfs);
    //! the YFS.NLO/NNLO breakdown columns; false if m_ev.m_nlo_current is off
    bool FillNLOBreakdownColumns(ATOOLS::Weights &wyfs, double w_lo, double w_full);
    bool BuildLadderWeights(ATOOLS::Weights &w);

    void ProbeSign(double w_full);
  };

}

#endif
