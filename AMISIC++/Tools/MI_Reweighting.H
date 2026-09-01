#ifndef AMISIC_Tools_MI_Reweighting_H
#define AMISIC_Tools_MI_Reweighting_H

#include "ATOOLS/Phys/Blob.H"
#include <vector>

namespace AMISIC {
  class MI_Parameters;
  class Matter_Overlap;
  class Interaction_Probability;
  class Hadronic_XSec_Calculator;
  class Over_Estimator;
  class MI_Processes;
  class Single_Collision_Handler;

  class MI_Reweighting {
  private:
    const MI_Parameters      * p_mipars;
    Matter_Overlap           * p_mo;
    Interaction_Probability  * p_pint;
    Hadronic_XSec_Calculator * p_xsecs;
    Over_Estimator           * p_overestimator;
    MI_Processes             * p_processes;
    Single_Collision_Handler * p_singlecollision;

    size_t m_n_variations;
    std::vector<double> m_sigma_nd_variations;
    std::vector<double> m_pt0_variations;
    std::vector<double> m_ptmin_variations;
    std::vector<double> m_eta_variations;
    std::vector<double> m_variation_weights;
    std::vector<double> m_b_weights;
    std::vector<double> m_overlap_ratios;
    std::vector<double> m_sudakov_weights;
    double m_max_reweight_factor;
  public:
    MI_Reweighting();
    ~MI_Reweighting();

    void Initialize(const MI_Parameters      * mipars,
                    Matter_Overlap           * mo,
                    Interaction_Probability  * pint,
                    Hadronic_XSec_Calculator * xsecs,
                    Over_Estimator           * overestimator,
                    MI_Processes             * processes,
                    Single_Collision_Handler * singlecollision);

    void ResetEvent();
    void SetOverlapRatios(const double & s, const double & b);
    void ImpactParameterReweighting(const double & s, const double & b,
                                    bool IsMinBias = false);
    void AcceptRejectReweighting(const bool accepted, const double prob_nom);
    void ApplyVariationWeights(ATOOLS::Blob * blob);
  };
}

#endif
