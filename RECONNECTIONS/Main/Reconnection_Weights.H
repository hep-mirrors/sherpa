#ifndef RECONNECTIONS_Main_Reconnection_Weights_H
#define RECONNECTIONS_Main_Reconnection_Weights_H

#include "ATOOLS/Phys/Particle.H"

namespace RECONNECTIONS {
  typedef std::map<ATOOLS::Particle *, double>  distances;

  class Reconnection_Base;
  
  class Reconnection_Weights {
  private:
    Reconnection_Base * p_reconnector;

    ATOOLS::ParticleSet                     * p_parts[2];
    std::map<ATOOLS::Particle *,distances * > m_distances;

    void   SetLists();
    void   OutputWeightTable();
  public:
    Reconnection_Weights(Reconnection_Base * reconnector);
    ~Reconnection_Weights() {}
    
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
