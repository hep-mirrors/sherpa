#ifndef RECONNECTIONS_Main_Reconnect_Statistical_H
#define RECONNECTIONS_Main_Reconnect_Statistical_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"

namespace RECONNECTIONS {
  class Reconnect_Statistical : public Reconnection_Base {
  private:
    std::vector<unsigned int> m_collist;
    int    m_Pmode;
    double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_kappa, m_norm;
    //double m_min_stringlength, m_max_stringlength; // temp ivars, will help inform appropriate hist. binning
    std::vector<double> m_stringlength_totals;

    void   SetParameters();
    bool   SelectColourPair(const size_t & N,unsigned int & col1, unsigned int & col2);
    bool   AttemptSwap(const unsigned int col[2]);
    void   UpdateColours();
    
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double TotalLength();
    void   PlotTotalLength(double total_length);
//    double find_min_totalLength(std::vector<double> total_vec);
//    double find_max_totalLength(std::vector<double> total_vec);
    void   write_stringLengths_to_file(std::vector<double> total_vec);
  public:
    Reconnect_Statistical();
    ~Reconnect_Statistical();

    int    operator()(ATOOLS::Blob_List *const blobs);
    void   Reset();
    double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  };
}
#endif
