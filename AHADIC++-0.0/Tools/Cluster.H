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


  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);

  struct ctp {
    enum code {
      no = 0,
      qq = 11,
      qd = 12,
      dq = 21,
      dd = 22
    };
  };

  class Cluster {
  private:
    ctp::code           m_type;
    ATOOLS::Flavour     m_flavours[2];
    ATOOLS::Vec4D       m_momentum, m_momenta[2];
    ATOOLS::Poincare    m_boost;
    bool                m_hasboost;

    Cluster           * p_left, * p_right;
  public:
    Cluster();
    Cluster(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,
	    const ATOOLS::Flavour &,const ATOOLS::Vec4D &);
    ~Cluster();


    Cluster         * GetLeft() const;
    Cluster         * GetRight() const;
    ATOOLS::Flavour   GetFlav(const int) const;

    void          SetLeft(Cluster *);
    void          SetRight(Cluster *);
    void          DeleteLeft();
    void          DeleteRight();
    void          Update();
    void          RescaleMomentum(ATOOLS::Vec4D);
    void          Boost(ATOOLS::Poincare &);
    void          BoostBack(ATOOLS::Poincare &);
    void          BoostInCMS();
    void          BoostBack();
    void          BoostBack(ATOOLS::Vec4D &);
    void          SetFlav(const int,const ATOOLS::Flavour &);
    void          SetMomentum(const int,const ATOOLS::Vec4D &);
    ATOOLS::Vec4D Momentum(const int=0) const;
    double        Mass2() const;
    double        Mass(const int=0) const;
    friend std::ostream& operator<<(std::ostream &,const Cluster &);
  };


  inline Cluster              * Cluster::GetLeft()  const { return p_left;  }
  inline Cluster              * Cluster::GetRight() const { return p_right; }
  inline void                   Cluster::SetLeft(Cluster * left)   { p_left  = left;  }
  inline void                   Cluster::SetRight(Cluster * right) { p_right = right; }
  inline void                   Cluster::DeleteLeft()     { delete p_left;  p_left  = NULL; }
  inline void                   Cluster::DeleteRight()    { delete p_right; p_right = NULL; }
  inline double                 Cluster::Mass2()    const { return m_momentum.Abs2(); }




  typedef std::list<Cluster*>          Cluster_List;
  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);
}
#endif
