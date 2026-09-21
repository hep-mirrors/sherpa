#ifndef RECONNECTIONS_Main_Reconnection_Reweighting_H
#define RECONNECTIONS_Main_Reconnection_Reweighting_H

#include "ATOOLS/Phys/Blob_List.H"
#include <vector>
#include <cstddef>

namespace RECONNECTIONS {
  class Reconnection_Reweighting {
  private:
    size_t m_n_variations;
    std::vector<double> m_etaQ2;
    std::vector<double> m_reshuffle;
    std::vector<double> m_variation_weights;
    double m_max_reweight_factor;
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
  };
}

#endif
