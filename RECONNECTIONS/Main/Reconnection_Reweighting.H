#ifndef RECONNECTIONS_Main_Reconnection_Reweighting_H
#define RECONNECTIONS_Main_Reconnection_Reweighting_H

#include "ATOOLS/Phys/Blob_List.H"
#include <vector>
#include <cstddef>
#include <fstream> // OUTPUT

namespace RECONNECTIONS {
  class Reconnection_Reweighting {
  private:
    size_t m_n_variations;
    std::vector<double> m_etaQ2;
    std::vector<double> m_reshuffle;
    std::vector<double> m_variation_weights;
    double m_max_reweight_factor;

    // OUTPUT
    bool m_reweighting_output;
    ATOOLS::Blob_List * p_blobs;
    std::vector<size_t> m_cutoff_count;
    std::vector<double> m_sum_weights;
    std::vector<double> m_sum_weights_squared;
    std::vector<double> m_event_variation_weights;
    size_t m_total_events;
    size_t m_n_reconnection_count;
    size_t m_event_n_reconnections;
    std::ofstream m_cr_weight_file;
    std::ofstream m_total_weight_file;
    void AccumulateEventStatistics();
    void PrintVariationStatistics();
  public:
    Reconnection_Reweighting();
    ~Reconnection_Reweighting();

    void Initialize();
    void ResetEvent();
    size_t NumVariations() const { return m_n_variations; }
    double EtaQ2(size_t i) const { return m_etaQ2[i]; }
    double Reshuffle(size_t i) const { return m_reshuffle[i]; }
    void AcceptRejectReweighting(bool accepted, const std::vector<double>& probs);
    void ApplyVariationWeights(ATOOLS::Blob_List *const blobs);

    // OUTPUT
    void ResetStats();
    void CacheBlobs(ATOOLS::Blob_List * blobs);
    ATOOLS::Blob_List * GetBlobs() const { return p_blobs; }
    void WriteEventStatistics();
  };
}

#endif
