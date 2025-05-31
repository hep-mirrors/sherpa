#ifndef NUCLEON_CLUSTERER_Main_Clusterer_Weights_H
#define NUCLEON_CLUSTERER_Main_Clusterer_Weights_H

#include "ATOOLS/Phys/Particle.H"

namespace NUCLEON_CLUSTERER {
  typedef std::map<ATOOLS::Particle *, double>  distances;

  class Clusterer_Base;
  
  class Clusterer_Weights {
  private:
    Clusterer_Base * p_clusterer;

    ATOOLS::ParticleSet                     * p_parts[2];
    std::map<ATOOLS::Particle *,distances * > m_distances;

    void   SetLists();
    void   OutputWeightTable();
  public:
    Clusterer_Weights(Clusterer_Base * reconnector);
    ~Clusterer_Weights() {}
    
    void FillTables();
    void Reset();

    inline const double operator()(ATOOLS::Particle * part1,
				   ATOOLS::Particle * part2) {
      if (part1==part2) return 1.;
      distances * dist = m_distances[part1];
      return (*dist)[part2];
    }
    inline void SetWeight(ATOOLS::Particle * part1,
			  ATOOLS::Particle * part2,const double & wt) {
      (*m_distances[part1])[part2] = wt;
    }
  };
}
#endif
