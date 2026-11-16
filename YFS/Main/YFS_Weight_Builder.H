#ifndef YFS_Main_YFS_Weight_Builder_H
#define YFS_Main_YFS_Weight_Builder_H

#include "ATOOLS/Phys/Weights.H"
#include <vector>
#include <string>
#include <utility>

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
      double corr_acraic_fo = 0., corr_acraic_fo2 = 0., corr_acraic_fo3 = 0.;
      //! the EEX correction of the event (N1_PROBE)
      double corr_eex = 0.;
      double corr_acraic_nob02 = 0., corr_acraic_novp = 0.;
      //! ACRAIC_RVNLL shift
      double rvnll_shift = 0.;
      double tnp_vhalf = 0., vext_fo = 0., dyif_fo = 0., dv2_acraic = 0.;
      //! TNP_IFI: the charge-odd part of v; FIXED_ORDER_STREAM ORDER 2:
      //! delta^(2) of BETA0_2 and the real-virtual remainder beyond (1 + v/2)
      double tnp_vif = 0., d2_fo = 0., rvrem_fo = 0.;
      double b02_acraic = 0., b02_v0 = 0., b02_v1 = 0.;
      double vp_acraic = 0., vp_v0 = 0., vp_v1 = 0.;
      std::vector<double> tnp1, tnp2, tnpifi;
      bool acraico2 = false, acraicb02 = false, acraicvp = false;
    };
    //! ACRAIC_Virtual: external/automatic - the squared-level provider's v on A_0
    void ApplyExternalVirtual(AcraicColumns &c);
    //! ACRAIC_Virtual: helicity - the same insertion, per helicity
    void ApplyHelicityVirtual(AcraicColumns &c);
    //! beta_0^(2) and Dyson VP, added after either virtual above; TNP templates
    void FinalizeAcraicBeta02VP(AcraicColumns &c);
    //! exp(Y_IF) cut at the generation cutoff (see GenerateWeight)
    void ApplyIFFormFactorCut(AcraicColumns &c);
    //! ACRAIC: FIXED_ORDER_STREAM
    void ComputeFixedOrderColumn(AcraicColumns &c);
    //! YFS: N1_PROBE - the weight pieces of a one-photon event (@@@ N1PROBE)
    void ProbeN1(const AcraicColumns &c, size_t nphot, double yfo);
    //! YFS: N0_PROBE: the zero-photon fixed-order soft + virtual pieces per event
    void ProbeN0(const AcraicColumns &c, size_t nphot, double yfo);
    //! ACRAIC: FIXED_ORDER_STREAM at ORDER 2: the O(alpha^2) truncation of
    //! exp(Y(emin)) rho on an event with nphot photons
    double AcraicFixedOrder2(size_t nphot, double y, double v, double d2,
                             double rvrem) const;
    //! ... one order up: plus the column's own O(alpha^3) Taylor coefficient
    double AcraicFixedOrder3(size_t nphot, double y, double v, double d2,
                             double rvrem) const;

    void AcraicCompare();
    void InitNLO();
    double CalculateNLO();
    double AcraicRealVirtualRemainder();
    //! ACRAIC: REAL_VIRTUAL isr_nll, the analytic LL+NLL ISR remainder
    double AcraicIsrNllRemainder();
    //! ACRAIC: REAL_VIRTUAL hybrid (isr_nll on ISR photons, the exact remainder on FSR photons)
    double AcraicHybridRemainder();
    double AcraicRealVirtualHelicityRemainder();
    //! ACRAIC's photon index -> YFS.NLO's real-virtual dv (and its m_rvdv entry);
    //! photons before first keep their dv
    int MatchRealVirtualDV(std::vector<double> &dv, size_t first = 0,
                           std::vector<long> *entry = nullptr);
    //! ACRAIC: ISR_VERTEX_AT_S: the initial vertex shift from the post-ISR Born
    //! point to s (the Born virtual), or to s_pt (an (n+1)-point's beams)
    double IsrVertexShift() const;
    double IsrVertexShift(double spt) const;
    //! ACRAIC: TNP_IFI: v_IF less greco's resonant log, which sits in the propagators
    double TnpVirtualIF() const;
    double AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const;
    //! ACRAIC: BETA0_2_ON beta1: delta^(2) on A_1 (Born-level and beta_1 terms)
    double AcraicBeta02RhoOnBeta1(double d2, double v) const;
    void BuildNamedWeights(double w_lo, double w_full);
    void RefuseStreamsWithoutBreakdown() const;
    //! the YFS.* column names this configuration emits, and wyfs/m_wnames defaults
    //! ACRAIC: TNP_IFI columns: name and index in m_ev.m_tnpifi
    typedef std::vector<std::pair<std::string, size_t> > TnpIfiColumns;
    void RegisterColumnNames(std::vector<std::string> &names,
                             std::vector<std::string> &tnpnames,
                             TnpIfiColumns &tnpifi,
                             ATOOLS::Weights &wyfs);
    void FillAcraicColumns(ATOOLS::Weights &wyfs);
    void FillTnpColumns(ATOOLS::Weights &wyfs, const std::vector<std::string> &tnpnames);
    void FillTnpIfiColumns(ATOOLS::Weights &wyfs, const TnpIfiColumns &tnpifi);
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
