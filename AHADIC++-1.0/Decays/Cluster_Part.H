#ifndef Cluster_Part_H
#define Cluster_Part_H

#include "Pair_Popper.H"
#include "Cluster.H"
#include "Histogram.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Cluster_Part {
  protected:
    int             m_4Qmode;
    bool            m_ana;
    ATOOLS::Vec4D   m_momenta[4];

    Pair_Popper * p_popper;
    double        m_fraction, m_ystar_expvalue, m_ystar_sigma;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool   BuildKinematics(Cluster * const,ATOOLS::Flavour &);
    double GetYStar(const double,const double);
  public:
    Cluster_Part();
    ~Cluster_Part();
    bool TestDecay(Cluster *);
    bool UpdateDecay(Cluster *,const int);
  };
}
#endif
