#include "AMISIC++/Tools/MI_Reweighting.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/Interaction_Probability.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>
#include <cmath>

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

MI_Reweighting::MI_Reweighting() :
  p_mipars(nullptr), p_mo(nullptr), p_pint(nullptr), p_xsecs(nullptr),
  p_overestimator(nullptr), p_processes(nullptr), p_singlecollision(nullptr),
  m_n_variations(1),
  m_max_reweight_factor(-1.)
{}

MI_Reweighting::~MI_Reweighting() {}

void MI_Reweighting::Initialize(const MI_Parameters      * mipars,
                                Matter_Overlap           * mo,
                                Interaction_Probability  * pint,
                                Hadronic_XSec_Calculator * xsecs,
                                Over_Estimator           * overestimator,
                                MI_Processes             * processes,
                                Single_Collision_Handler * singlecollision) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialize the parameter variations for MPI reweighting.
  ///////////////////////////////////////////////////////////////////////////
  p_mipars          = mipars;
  p_mo              = mo;
  p_pint            = pint;
  p_xsecs           = xsecs;
  p_overestimator   = overestimator;
  p_processes       = processes;
  p_singlecollision = singlecollision;

  if (!p_mipars) return;
  m_sigma_nd_variations = p_mipars->GetVariationVector("SigmaND_Norm");
  m_pt0_variations      = p_mipars->GetVariationVector("pt_0");
  m_ptmin_variations    = p_mipars->GetVariationVector("pt_min");
  m_eta_variations      = p_mipars->GetVariationVector("eta");
  size_t n_variations_without_matterform = std::max({m_sigma_nd_variations.size(),
                                                     m_pt0_variations.size(),
                                                     m_ptmin_variations.size(),
                                                     m_eta_variations.size(),
                                                     size_t(1)});
  m_n_variations = std::max(n_variations_without_matterform,
                              p_mo->MatterFormVariationSize());

  m_pt0_variations.resize(m_n_variations, m_pt0_variations[0]);
  m_ptmin_variations.resize(m_n_variations, m_ptmin_variations[0]);
  if (m_n_variations > 1) {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      if (m_ptmin_variations[ivar]<m_ptmin_variations[0]) {
        THROW(fatal_error, std::string("Reweighting of MPI only possible for upward variations of PT_Min.\n")
                          + "Found PT_Min variation " + std::to_string(ivar) + " = "
                          + std::to_string(m_ptmin_variations[ivar]) + " < "
                          + std::to_string(m_ptmin_variations[0]) + " = PT_Min nominal.\n"
                          + "Please adjust your parameter variation settings (PT_Min or/and Eta).");
      }
    }
  }

  m_sigma_nd_variations.resize(n_variations_without_matterform, m_sigma_nd_variations[0]);

  // Set relevant variations for total cross section computation (matter form independent)
  p_xsecs->SetVariations(n_variations_without_matterform, m_sigma_nd_variations);

  if (m_n_variations > n_variations_without_matterform) {
  m_sigma_nd_variations.resize(m_n_variations, m_sigma_nd_variations[0]); }

  // Set relevant variations for interaction probability and overestimator (matter form dependent)
  p_pint->SetVariations(m_n_variations, m_sigma_nd_variations);
  p_overestimator->SetVariations(m_n_variations);

  m_max_reweight_factor = (*p_mipars)("max_reweight_factor");
  ResetEvent();

  if (m_n_variations > 1)
    RegisterSoftPhysicsVariations(m_n_variations);
}

void MI_Reweighting::CheckReweightingSampling(const double & s, bool IsMinBias) {
  ///////////////////////////////////////////////////////////////////////////
  // Impact parameters are only drawn from the nominal density: warn for any
  // variation with more than 1% of its b density beyond the nominal 99.9%
  // quantile, where the nominal sample has (almost) no support.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1 || !p_pint || !p_mo) return;
  const axis bbins(4000, 0., p_pint->Bmax(), axis_mode::linear);
  double norm_nom = 0., norm_var = 0.;
  const std::unique_ptr<OneDim_Table> cnom =
    BSamplingDensity(s, bbins, 0, IsMinBias).Cumulative(1., norm_nom);
  if (!(norm_nom>0.)) return;
  const double b999 = cnom->Inverse(0.999*norm_nom);
  const double knom = p_pint->K(s);
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    const std::unique_ptr<OneDim_Table> cvar =
      BSamplingDensity(s, bbins, ivar, IsMinBias).Cumulative(1., norm_var);
    const double cover = (norm_var>0.) ? (*cvar)(b999)/norm_var : 0.;
    if (cover < 0.99)
      msg_Error()<<"Warning in "<<METHOD<<": nominal b sampling covers only "
                 <<(100.*cover)<<"% of variation v"<<ivar<<" (k_var/k_nom = "
                 <<(knom>0. ? p_pint->K(s,ivar)/knom : 0.)
                 <<"), expect large weights in the peripheral b tail.\n";
  }
}

OneDim_Table MI_Reweighting::BSamplingDensity(const double & s,
                                              const axis & bbins,
                                              const size_t ivar,
                                              const bool IsMinBias) {
  ///////////////////////////////////////////////////////////////////////////
  // Tabulate the density per unit area b is sampled from for variation ivar:
  // P_int(s,b) for MinBias and the matter overlap O(b) otherwise.
  ///////////////////////////////////////////////////////////////////////////
  OneDim_Table density(bbins);
  if (!IsMinBias) p_mo->SetMatterFormVariationIndex(ivar);
  for (size_t i=0; i<=bbins.m_nbins; ++i) {
    const double b = bbins.x(i);
    const double q = ( IsMinBias ?
                       (ivar==0 ? (*p_pint)(s,b) : (*p_pint)(s,b,ivar)) :
                       (*p_mo)(b, ivar==0 ? p_pint->K(s) : p_pint->K(s,ivar)) );
    density.Fill(i, (std::isfinite(q) && q>0.) ? q : 0.);
  }
  p_mo->SetMatterFormVariationIndex(0);
  return density;
}

void MI_Reweighting::InitializeMinBiasNormalisation(const axis * sbins) {
  ///////////////////////////////////////////////////////////////////////////
  // Correct MinBias normalisation for discarded impact-parameter attempts.
  // Each b is retained with probability A = 1 - (1-P_int)^N, where N is the
  // number of scatter tries. The correction is the ratio of nominal to
  // varied integrals of P_int*A, each normalised by its sigma_ND.
  // Applied per energy and variation after P_int initialisation.
  ///////////////////////////////////////////////////////////////////////////
  m_minbias_norm.clear();
  if (m_n_variations <= 1 || !p_pint || !sbins) return;
  MinBias_Kept_Integrand integrand(p_pint,Single_Collision_Handler::s_minbias_tries);
  Gauss_Integrator       gauss(&integrand);

  for (size_t ivar=0; ivar<m_n_variations; ++ivar) m_minbias_norm.emplace_back(*sbins);
  std::vector<double> kept(m_n_variations, 0.);
  const size_t nnodes = sbins->m_nbins > 1 ? sbins->m_nbins : 0;
  for (size_t sbin=0; sbin<=nnodes; ++sbin) {
    const double s = sbins->x(sbin);
    integrand.SetS(s);
    for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
      integrand.SetVariation(ivar);
      kept[ivar] = s<=4. ? 0. : gauss.Integrate(0.,p_pint->Bmax(),1.e-6);
    }
    for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
      double c = ( (kept[0]/m_sigma_nd_variations[0]) /
                   (kept[ivar]/m_sigma_nd_variations[ivar]) );
      if (!std::isfinite(c) || c<=0.) c = 1.;
      m_minbias_norm[ivar].Fill(sbin, c);
    }
  }
}

double MI_Reweighting::MinBiasNormalisation(const double & s,
                                            const size_t ivar) const {
  if (ivar>=m_minbias_norm.size()) return 1.;
  const axis & sbins = m_minbias_norm[ivar].GetAxis();
  const double c = m_minbias_norm[ivar](Max(sbins.m_xmin,Min(s,sbins.m_xmax)));
  return (std::isfinite(c) && c>0.) ? c : 1.;
}

double MinBias_Kept_Integrand::operator()(double b) {
  ///////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b P_int(b) A(b) = 2 pi b db P_int(b) A(b), where
  // A(b) = 1 - (1-P_int(b))^N is the probability that b is retained in
  // MinBias events, i.e. not discarded after N failed scatter tries.
  ///////////////////////////////////////////////////////////////////////////
  const double pint = m_ivar==0 ? (*p_pint)(m_s,b) : (*p_pint)(m_s,b,m_ivar);
  return 2.*M_PI*b*pint*(1.-exp(m_ntries*log1p(-pint)));
}

void MI_Reweighting::ResetEvent() {
  ///////////////////////////////////////////////////////////////////////////
  // Reset variation weights to 1.0 for a new event.
  ///////////////////////////////////////////////////////////////////////////
  m_variation_weights.resize(m_n_variations);
  std::fill(m_variation_weights.begin(), m_variation_weights.end(), 1.);
  m_bsampling_weights.assign(m_n_variations, 1.);
  m_overlap_ratios.assign(m_n_variations, 1.);
  m_sudakov_weights.assign(m_n_variations, 1.);
}

void MI_Reweighting::SetOverlapRatios(const double & s, const double & b) {
  ///////////////////////////////////////////////////////////////////////////
  // Set m_overlap_ratios[ivar] = O_var(b) / O_nom(b) at the (s, b) point.
  // Used by AcceptRejectReweighting for the Sudakov probability ratio, and
  // by ImpactParameterReweighting. Called standalone from FirstRescatter,
  // which inherits b from the prior event and must NOT re-apply a b-weight,
  // and from SelectPT2 for dynamic overlaps, whose ratio changes with the
  // x of every trial.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  const double overlap_nom = (*p_mo)(b);
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    p_mo->SetMatterFormVariationIndex(ivar);
    const double overlap_var = (*p_mo)(b, p_pint->K(s, ivar));
    double overlap_ratio = overlap_var / overlap_nom;
    if (!std::isfinite(overlap_ratio) || overlap_ratio<=0.) overlap_ratio = 1.;
    m_overlap_ratios[ivar] = overlap_ratio;
  }
  p_mo->SetMatterFormVariationIndex(0);
}

void MI_Reweighting::ImpactParameterReweighting(const double & s, const double & b,
                                                bool IsMinBias) {
  ///////////////////////////////////////////////////////////////////////////
  // Cache per-event information needed for Sudakov reweighting (via
  // SetOverlapRatios) and compute the impact-parameter weight w_b for the
  // b-sampling distribution of this event type:
  //   MPI     : b ~ d^2b O(b)
  //             w_b *= O_var(b) / O_nom(b)
  //   MinBias : b ~ d^2b P_int/sigma_ND
  //             w_b  = [P_int_var/sigma_ND_var] / [P_int_nom/sigma_ND_nom]
  //                    * C_var(s), see InitializeMinBiasNormalisation
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  SetOverlapRatios(s, b);
  if (IsMinBias) {
    const double pint_nom = (*p_pint)(s, b);
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      double w_b = 1.;
      const double pint_var = (*p_pint)(s, b, ivar);
      double pint_ratio = pint_var / pint_nom;
      w_b = ( pint_ratio * m_sigma_nd_variations[0] / m_sigma_nd_variations[ivar]
              * MinBiasNormalisation(s, ivar) );
      if (!std::isfinite(w_b) || w_b < 0.) w_b = 1.;
      m_bsampling_weights[ivar] *= w_b;
    }
  } else {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      m_bsampling_weights[ivar] *= m_overlap_ratios[ivar];
    }
  }
}

void MI_Reweighting::AcceptRejectReweighting(const bool accepted, const double prob_nom) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Single_Collision_Handler for each accept/reject event
  // during Sudakov evolution of MPI scatters.
  // For accepted scatter: multiply weight by p_var/p_nom
  // For rejected scatter: multiply weight by (1-p_var)/(1-p_nom)
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  if (prob_nom <= 1e-8) return;

  // PT_0 enters dsigma/dpt^2 only through Coupling = alpha_s^2(Max(pt0^2, m_muR_fac*(pt^2+pt0^2)))
  // and SoftCorrection = (pt^2/(pt^2+pt0^2))^2.
  const double pt2      = p_singlecollision->PT2();
  const double pt02_nom = sqr(m_pt0_variations[0]);
  MODEL::One_Running_AlphaS * alphaS = p_processes->AlphaS();
  const double muR_fac = p_processes->MuRFac();

  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    // p_var/p_nom = overlap_ratio * xs_ratio
    // p_var = 0 for pt2 <= (ptmin2)_var
    double prob_var = 0.;
    if (pt2 > sqr(m_ptmin_variations[ivar])) {
      double xs_ratio = 1.;
      if (m_pt0_variations[ivar] != m_pt0_variations[0]) {
        const double pt02_var  = sqr(m_pt0_variations[ivar]);
        const double scale_nom = Max(pt02_nom, muR_fac * (pt2 + pt02_nom));
        const double scale_var = Max(pt02_var, muR_fac * (pt2 + pt02_var));
        const double alpha_nom = (*alphaS)(scale_nom);
        const double alpha_var = (*alphaS)(scale_var);
        xs_ratio = sqr(alpha_var / alpha_nom)
                 * sqr((pt2 + pt02_nom) / (pt2 + pt02_var));
      }
      prob_var = prob_nom * m_overlap_ratios[ivar] * xs_ratio;
    }
    if (accepted) {
      // Accepted: weight *= p_var / p_nom
      const double ratio = prob_var / prob_nom;
      if (std::isfinite(ratio) && ratio >= 0.) {
        m_sudakov_weights[ivar] *= ratio;
      }
    } else {
      // Rejected: weight *= (1 - p_var) / (1 - p_nom)
      const double ratio = (1. - prob_var) / (1. - prob_nom);
      if (std::isfinite(ratio) && ratio >= 0.) {
        m_sudakov_weights[ivar] *= ratio;
      }
    }
  }
}

void MI_Reweighting::ApplyVariationWeights(ATOOLS::Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute and apply variation weights to the event.
  // The total weight is:
  // w_total = w_b * w_(n|b) = w_b * w_sudakov
  // where w_sudakov accounts for the changed MPI multiplicity distribution
  // and w_b accounts for the changed impact parameter sampling distribution.
  ///////////////////////////////////////////////////////////////////////////
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    double w_total = m_bsampling_weights[ivar] * m_sudakov_weights[ivar];
    if (m_max_reweight_factor > 0. && w_total > m_max_reweight_factor) {
      w_total = m_max_reweight_factor;
    }
    m_variation_weights[ivar] = w_total;
  }
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    CombineSoftPhysicsVariations(wgtmap, m_variation_weights);
    blob->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
  }
  ResetEvent();
}
