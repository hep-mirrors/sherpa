#ifndef Cluster_Part_H
#define Cluster_Part_H

#include "Pair_Popper.H"
#include "Cluster.H"


namespace AHADIC {
  class Cluster_Part {
  public:
    Cluster_Part() {}
    ~Cluster_Part() {}
    virtual bool TestDecay(Cluster *,Part_List *)=NULL;
  };

  class QoverM : public Cluster_Part {
    Pair_Popper m_popper;
    double      m_Q;
  public:
    QoverM();
    ~QoverM();
    bool TestDecay(Cluster *,Part_List *);
  };

  class FourFermion : public Cluster_Part {
  public:
    FourFermion() {}
    ~FourFermion() {}
    bool TestDecay(Cluster *,Part_List *) {}
  };
}
#endif
