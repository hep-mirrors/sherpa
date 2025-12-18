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

    void   SetParameters();
    void   FixPMode(const std::string & pm);
    bool   PerformTableReconnections(const size_t & N);
    void   UpdateColours();
    
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
  };
}
#endif
