#ifndef AHADIC_Tools_Cluster_H
#define AHADIC_Tools_Cluster_H

#include "AHADIC++/Tools/Proto_Particle.H"
#include <list>
#include <iostream>

namespace AHADIC {
  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);
  
  class Cluster {
  private:
    std::pair<Proto_Particle *,Proto_Particle *> m_parts;
    ATOOLS::Vec4D m_momentum;
  public:
    Cluster(std::pair<Proto_Particle *,Proto_Particle *> parts);
    Cluster(Proto_Particle * part1,Proto_Particle * part2);
    ~Cluster();
    void Clear();
    ATOOLS::Vec4D Momentum() { return m_momentum; } 
    Proto_Particle * operator[](size_t pos) {
      if (pos==0) return m_parts.first;
      if (pos==1) return m_parts.second;
      return NULL;
    }
    
    friend std::ostream & operator<<(std::ostream &,const Cluster &);
  };
  
  typedef std::list<Cluster *>         Cluster_List;
  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;
  
  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);
}
#endif
