#ifndef AHADIC_Tools_Hadronisation_Reweighting_H
#define AHADIC_Tools_Hadronisation_Reweighting_H

#include "ATOOLS/Phys/Blob.H"
#include <vector>
#include <map>
#include <string>
// OUTPUT
#include "ATOOLS/Phys/Flavour.H"
#include <utility>
#include <fstream>

namespace ATOOLS { class Histogram; }

namespace AHADIC {

  class Hadronisation_Reweighting {
  public:

    /*!
      Logarithm of the normalisation of the cluster-splitting fragmentation
      function,

        N(alpha,beta,c) = int_zmin^zmax dz z^alpha (1-z)^beta exp(-c/z),

      as needed by the on-the-fly reweighting.  Only the reweighting uses this;
      the z sampling in Cluster_Splitter::FragmentationFunction never does, so
      nothing in here can change the generated event sample or the random
      number sequence.

      The node set is built once per accepted z, by SetRange(), and is then
      common to every variation.  That is deliberate: the reweighting needs
      N_var/N_nom, and on a common node set the quadrature errors of the two
      integrals largely cancel.  For the same reason the rule must not be
      replaced by one whose nodes depend on alpha, beta or c.

      Resolution.  The exponent alpha*ln z + beta*ln(1-z) - c/z is linear in
      u = 1/z, and the exp(-c/z) factor confines the integrand to a window of
      width ~1/c in u around its peak.  A grid graded only in z and in 1-z has
      panel widths that do not depend on c at all, so beyond some c the peak
      falls between two nodes and N is underestimated by an unbounded factor.
      Because the numerator of the reweighting factor is evaluated exactly,
      that error goes straight into the weight, and it does not cancel between
      variations with different c.  The grid is therefore the union of

        - a base grid graded geometrically in z below the centre of the
          interval and in 1-z above it, which resolves z^alpha, (1-z)^beta and
          any moderate c, and

        - when the base grid is too coarse for the largest c in play, a
          refinement graded geometrically in |u - u_peak| about the peak of
          the most strongly peaked variation.  Its first step resolves the
          largest c; the geometric growth then resolves each smaller c out at
          the distance where that variation is still significant.

      The refinement is skipped entirely when the base grid already resolves
      the largest c, which is the common case, so ordinary splittings pay
      nothing.
    */
    class Frag_Norm {
    public:
      Frag_Norm();

      /*!
        Prepare the node set for [zmin,zmax].  alpha and beta are those of the
        most strongly peaked variation and cmin/cmax the smallest and largest
        |c| over all variations; they only steer where nodes are placed, never
        what is integrated.  Requires 0 < zmin < zmax < 1.
      */
      void SetRange(const double zmin, const double zmax,
                    const double alpha, const double beta,
                    const double cmin, const double cmax);

      //! log N(alpha,beta,c) on the node set prepared by SetRange().
      double operator()(const double alpha, const double beta,
                        const double c) const;

      size_t NNodes()  const { return m_nnodes; }
      size_t NPanels() const { return m_edge.empty() ? 0 : m_edge.size()-1; }

      //! Peak of z^alpha (1-z)^beta exp(-c/z) on [zmin,zmax].
      static double Peak(const double alpha, const double beta, const double c,
                         const double zmin, const double zmax);

      /*!
        The exponent c of the fragmentation function, shared with
        Cluster_Splitter::FragExponent so that the self-check below and the
        reweighting cannot drift apart.
      */
      static double Exponent(const double gamma, const double kt02,
                             const double scale);

      /*!
        Startup self-check, run once from Initialize() for the variations that
        are actually configured.  Sweeps the kinematic reach of the cluster
        splitting and compares this rule against a deliberately much finer one.
        Silent unless it finds a problem.
      */
      static void Verify(const Hadronisation_Reweighting& reweighting);

      /*!
        Grid controls.  The defaults are the ones used in production; Verify()
        builds a second, deliberately finer instance and compares against it.
        maxdphi  - exponent change tolerated across one panel
        grading  - growth factor of the refinement steps
        refrange - how far down from the peak, in units of the exponent, the
                   refinement is carried
        panels   - hard cap on the number of panels
      */
      void SetGrid(const double maxdphi, const double grading,
                   const double refrange, const int panels);

    private:
      std::vector<double> m_w, m_logz, m_log1mz, m_invz, m_edge;
      mutable std::vector<double> m_e;
      size_t m_nnodes;
      double m_maxdphi, m_grading, m_refrange;
      int    m_maxpanels;

      void BuildNodes();
    };

  private:
    size_t m_n_variations;
    size_t m_n_flavour_variations;
    size_t m_n_gluon_variations;
    size_t m_n_cluster_variations;
    size_t m_n_soft_variations;
    size_t m_n_kt_variations;
    size_t m_n_promptdecay_variations;
    double m_max_reweight_factor;
    int    m_reweight_max_nsplit;
    bool   m_frag_ar {false};
    std::map<std::string, std::vector<double>> m_variation_vectors;

    std::vector<double> m_variation_weights;
    // Per-channel event weights, accumulated over all splittings of an event:
    // flavour popping (Flavour_Selector), gluon splitting z (Gluon_Splitter),
    // cluster splitting z (Cluster_Splitter), soft cluster decays and
    // prompt-decay decisions (Soft_Cluster_Handler), and transverse
    // momentum (KT_Selector).
    std::vector<double> m_flavour_weights;
    std::vector<double> m_gluon_weights;
    std::vector<double> m_cluster_weights;
    std::vector<double> m_soft_weights;
    std::vector<double> m_kt_weights;
    std::vector<double> m_promptdecay_weights;
    // Per-splitting weights, accumulated over every attempt of a splitting
    // and committed into the event weights when the splitting succeeds.
    std::vector<double> m_tmp_flavour_weights;
    std::vector<double> m_tmp_gluon_weights;
    std::vector<double> m_tmp_cluster_weights;
    // The remaining three channels are committed directly while no splitting
    // is open, and buffered like the ones above while one is.
    bool m_in_splitting {false};
    std::vector<double> m_tmp_kt_weights;
    std::vector<double> m_tmp_soft_weights;
    std::vector<double> m_tmp_promptdecay_weights;

    void CheckVariationGuards();
    void AcceptTmpWeights(size_t n_variations, std::vector<double>& weights,
                          std::vector<double>& tmp_weights);
    void ResetTmpWeights(size_t n_variations, std::vector<double>& tmp_weights);

    // OUTPUT
    int  m_output_mode;
    bool m_reweighting_output;
    bool m_event_weights_applied;
    size_t m_total_events;
    std::vector<double> m_event_variation_weights;
    std::vector<double> m_event_flavour_weights;
    std::vector<double> m_event_gluon_weights;
    std::vector<double> m_event_cluster_weights;
    std::vector<double> m_event_soft_weights;
    std::vector<double> m_event_kt_weights;
    std::vector<double> m_event_promptdecay_weights;
    std::vector<double> m_sum_weights;
    std::vector<double> m_sum_weights_squared;
    std::vector<size_t> m_cutoff_count;
    void AccumulateEventStatistics();
    void PrintVariationStatistics();
    std::ofstream m_hadronisation_weight_file;
    typedef std::vector<std::pair<char, std::string> > Record_Vector;
    Record_Vector m_tmp_flavour_records;
    Record_Vector m_tmp_gluon_records;
    Record_Vector m_tmp_cluster_records;
    Record_Vector m_call_records;
    Record_Vector m_event_records;
    void CommitTmpRecords(Record_Vector& tmp_records);
    void DiscardTmpRecords(Record_Vector& tmp_records);
    std::string m_hist_path;
    std::map<std::string, ATOOLS::Histogram*> m_hist_nominal;
    std::map<std::string, std::vector<ATOOLS::Histogram*> > m_hist_reweighted;
    void BookHistograms();
    void FillHistograms();
    void FillWeightHistogram();
    void WriteHistograms();
    void FillObs(const std::string& key, const double value);
  public:
    Hadronisation_Reweighting();
    ~Hadronisation_Reweighting();

    void Initialize();
    void ResetEvent();

    inline bool   Active() const { return m_n_variations > 1; }
    inline bool   OutputActive() const { return m_reweighting_output > 0; }
    inline size_t NumberOfVariations() const { return m_n_variations; }
    inline size_t NumberOfFlavourVariations() const { return m_n_flavour_variations; }
    inline size_t NumberOfGluonVariations() const { return m_n_gluon_variations; }
    inline size_t NumberOfClusterVariations() const { return m_n_cluster_variations; }
    inline size_t NumberOfSoftVariations() const { return m_n_soft_variations; }
    inline size_t NumberOfKTVariations() const { return m_n_kt_variations; }
    inline size_t NumberOfPromptDecayVariations() const { return m_n_promptdecay_variations; }
    std::vector<double> GetVariationVector(const std::string& keyword) const;

    void CheckTransitionGuard(const ATOOLS::Flavour& hadron,
                              const std::vector<double>& weights,
                              const double& threshold) const;
    void CheckPoppingGuard(const ATOOLS::Flavour& popped,
                           const std::vector<double>& weights,
                           const double& threshold) const;
    void CheckComponentGuard(const ATOOLS::Flavour& hadron,
                             const std::vector<double>& amps,
                             const double& threshold) const;

    // Flavour popping in the shared Flavour_Selector
    void FlavourSelectionReweighting(const std::vector<double>& popweights,
                                     const std::vector<double>& norms);
    inline void AcceptFlavourSelectionWeights() {
      AcceptTmpWeights(m_n_flavour_variations, m_flavour_weights, m_tmp_flavour_weights);
      CommitTmpRecords(m_tmp_flavour_records); // OUTPUT
    }
    inline void ResetFlavourSelectionWeights() {
      ResetTmpWeights(m_n_flavour_variations, m_tmp_flavour_weights);
      DiscardTmpRecords(m_tmp_flavour_records); // OUTPUT
    }

    // Gluon splitting z selection (analytic density ratio per drawn z)
    void GluonSplittingReweighting(const std::vector<double>& probs);
    inline void AcceptGluonSplittingWeights() {
      AcceptTmpWeights(m_n_gluon_variations, m_gluon_weights, m_tmp_gluon_weights);
      CommitTmpRecords(m_tmp_gluon_records); // OUTPUT
    }
    inline void ResetGluonSplittingWeights() {
      ResetTmpWeights(m_n_gluon_variations, m_tmp_gluon_weights);
    }
    inline void DiscardGluonSplittingRecords() { // OUTPUT
      DiscardTmpRecords(m_tmp_gluon_records);
    }

    // Cluster splitting z selection (analytic density ratio per drawn z)
    inline bool DoClusterSplittingReweighting(const int nsplit) const {
      return m_n_cluster_variations > 1 &&
             (m_reweight_max_nsplit < 0 || nsplit < m_reweight_max_nsplit);
    }
    inline bool ARFragReweighting() const { return m_frag_ar; }
    void ClusterSplittingReweighting(const std::vector<double>& logprobs);
    void ClusterSplittingReweightingAR(const bool accepted,
                                       const std::vector<double>& probs);
    inline void AcceptClusterSplittingWeights() {
      AcceptTmpWeights(m_n_cluster_variations, m_cluster_weights, m_tmp_cluster_weights);
      CommitTmpRecords(m_tmp_cluster_records); // OUTPUT
    }
    inline void ResetClusterSplittingWeights() {
      ResetTmpWeights(m_n_cluster_variations, m_tmp_cluster_weights);
    }
    inline void DiscardClusterSplittingRecords() { // OUTPUT
      DiscardTmpRecords(m_tmp_cluster_records);
    }

    void BeginSplitting();
    void CommitSplitting();
    void AbortSplitting();

    // Transverse momentum selection
    void KTSelectionReweighting(const std::vector<double>& probs);

    // Hadron pair selection in soft cluster decays
    void SoftClusterReweighting(const std::vector<double>& weights,
                                const std::vector<double>& totweights);

    // Bernoulli prompt-decay decision, only stochastic for zeta > 0
    void PromptDecayReweighting(const bool decayed,
                                const std::vector<double>& probs);

    void ApplyVariationWeights(ATOOLS::Blob * blob);

    // OUTPUT
    void ResetStats();
    void WriteEventStatistics();
    void RecordFlavourPop(const ATOOLS::Flavour& flav, const double& Emax);
    void RecordGluonZ(const double& z);
    void RecordClusterZ(const int nsplit,
                        const int type1, const double& z1,
                        const int type2, const double& z2,
                        const double& mass);
    void RecordKT(const int site, const double& kt, const double& ktmax);
    void RecordSoftDecay(const double& mass, const ATOOLS::Flavour& had1,
                         const ATOOLS::Flavour& had2);
    void RecordTransition(const double& mass, const ATOOLS::Flavour& had);
    void RecordPrimaryCluster(const double& mass, const ATOOLS::Flavour& fl1,
                              const ATOOLS::Flavour& fl2, const bool direct);
  };
}

#endif
