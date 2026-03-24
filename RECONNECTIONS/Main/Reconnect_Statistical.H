#ifndef RECONNECTIONS_Main_Reconnect_Statistical_H
#define RECONNECTIONS_Main_Reconnect_Statistical_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"

namespace RECONNECTIONS {
  class Reconnect_Statistical : public Reconnection_Base {
  private:
    std::vector<unsigned int> m_collist;
    int    m_Pmode;
    double m_Q02, m_R02, m_etaR, m_kappa, m_norm;
    std::vector<double> m_etaQ, m_reshuffle;
    double m_max_reweight_factor;
    bool m_weight_output;

    void   SetParameters();
    void   FixPMode(const std::string & pm);
    bool   PerformTableReconnections(const size_t & N);
    void   UpdateColours();
    size_t m_n_variations;

    size_t m_n_reconnection_count;
    double m_initial_length;
    double m_final_length;
    
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double TotalLength();
  public:
    Reconnect_Statistical();
    ~Reconnect_Statistical();

    int    operator()(ATOOLS::Blob_List *const blobs);
    void   Reset();
    double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double GetMaxReweightFactor() const { return m_max_reweight_factor; }
    bool   GetWeightOutput() const { return m_weight_output; }
    int    GetVariationSize() const { return m_n_variations; }
    size_t GetReconnectionCount() const { return m_n_reconnection_count; }
    double GetInitialLength() const { return m_initial_length; }
    double GetFinalLength() const { return m_final_length; }
  };
}
#endif
