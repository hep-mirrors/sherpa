#ifndef Cluster_H
#define Cluster_H

#include "Particle.H"
#include "Poincare.H"
#include "Message.H"
#include <list>
#include <iostream>

namespace AHADIC {
  typedef std::list<ATOOLS::Particle*> Part_List;
  typedef Part_List::iterator Part_Iterator;
  typedef Part_List::const_iterator Part_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Part_List & pl);

  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour> FlavourPair;
  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);


  class Cluster {
  private:
    int                m_type;
    FlavourPair      * p_fpair;
    ATOOLS::Particle * p_parts[2];
    ATOOLS::Vec4D      m_momentum;
    ATOOLS::Poincare   m_boost;
    bool               m_hasboost;
  public:
    Cluster();
    Cluster(ATOOLS::Particle * trip, ATOOLS::Particle * antitrip);
    ~Cluster();

    Cluster          * p_left, * p_right;

    FlavourPair      * FlavPair();
    ATOOLS::Particle * GetParticle(const int i);

    void          RescaleMomentum(ATOOLS::Vec4D);
    void          Boost(ATOOLS::Poincare &);
    void          BoostBack(ATOOLS::Poincare &);
    void          BoostInCMS();
    void          BoostBack();
    void          BoostBack(ATOOLS::Vec4D &);
    void          SetFlavour(ATOOLS::Flavour flav);
    void          SetMomentum(ATOOLS::Vec4D mom);
    ATOOLS::Vec4D Momentum(const int i=-1);
    double        Mass2();
    double        Mass(const int i=-1);
    friend std::ostream& operator<<(std::ostream &,const Cluster &);
  };

  inline FlavourPair * Cluster::FlavPair() { return p_fpair; }
  inline double Cluster::Mass2() { return m_momentum.Abs2(); }
  inline void Cluster::BoostBack(ATOOLS::Vec4D & mom) {
    if (!m_hasboost) return;
    m_boost.BoostBack(mom);
  }

  inline void               Cluster::SetFlavour(ATOOLS::Flavour flav) { 
    p_fpair->first = flav; p_fpair->second = ATOOLS::Flavour(ATOOLS::kf::none); 
  }
  inline void               Cluster::SetMomentum(ATOOLS::Vec4D mom) { 
    m_momentum = mom;
  }



  typedef std::list<Cluster*>          Cluster_List;
  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);
}
#endif
