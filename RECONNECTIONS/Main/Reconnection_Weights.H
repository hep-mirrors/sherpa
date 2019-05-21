#ifndef RECONNECTIONS_Main_Reconnection_Weights_H
#define RECONNECTIONS_Main_Reconnection_Weights_H

#include "ATOOLS/Phys/Particle.H"

namespace RECONNECTIONS {
  typedef std::map<ATOOLS::Particle *, double>  distances;

  class Reconnection_Handler;
  
  class Reconnection_Weights {
  private:
    Reconnection_Handler * p_rhandler;
    int    m_Pmode;
    double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_restring;

    ATOOLS::ParticleSet            * p_parts[2];
    std::list<ATOOLS::Part_List *> * p_singlets;
    std::map<ATOOLS::Particle *,distances * > m_distances;

    void   SetLists();
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    void   OutputWeightTable();
  public:
    Reconnection_Weights(Reconnection_Handler * rhandler);
    ~Reconnection_Weights() {}

    void Initialize();
    void FillTables();
    void Reset();

    inline const double operator()(ATOOLS::Particle * part1,ATOOLS::Particle * part2) {
      if (part1==part2) return 1.;
      distances * dist = m_distances[part1];
      return (*dist)[part2];
    }
    inline void SetWeight(ATOOLS::Particle * part1,ATOOLS::Particle * part2,const double & wt) {
      (*m_distances[part1])[part2] = wt;
    }
  };
}
#endif
