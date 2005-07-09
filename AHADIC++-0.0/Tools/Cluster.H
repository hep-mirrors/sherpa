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

    Cluster           * p_left, * p_right, * p_prev;
  public:
    Cluster();
    Cluster(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,
	    const ATOOLS::Flavour &,const ATOOLS::Vec4D &);
    ~Cluster();


    Cluster         * GetLeft()  const;
    Cluster         * GetRight() const;
    Cluster         * GetPrev()  const;
    void              SetLeft(Cluster *);
    void              SetRight(Cluster *);
    void              SetPrev(Cluster *);
    void              DeleteLeft();
    void              DeleteRight();
    ATOOLS::Flavour   GetFlav(const int) const;
    void              SetFlav(const int,const ATOOLS::Flavour &);
    ATOOLS::Vec4D     Momentum(const int=0) const;
    void              SetMomentum(const int,const ATOOLS::Vec4D &);
    double            Mass2() const;
    double            Mass(const int=0) const;
    void              Update();
    void              RescaleMomentum(ATOOLS::Vec4D);
    void              BoostInCMS();
    void              BoostBack();
    void              Boost(ATOOLS::Poincare &);
    void              BoostBack(ATOOLS::Poincare &);
    void              BoostBack(ATOOLS::Vec4D &);
    friend std::ostream& operator<<(std::ostream &,const Cluster &);
  };


  inline Cluster              * Cluster::GetLeft()  const { return p_left;  }
  inline Cluster              * Cluster::GetRight() const { return p_right; }
  inline Cluster              * Cluster::GetPrev()  const { return p_prev;  }
  inline void                   Cluster::SetLeft(Cluster * left)   { p_left  = left;  }
  inline void                   Cluster::SetRight(Cluster * right) { p_right = right; }
  inline void                   Cluster::SetPrev(Cluster * prev)   { p_prev  = prev;  }
  inline void                   Cluster::DeleteLeft()     { delete p_left;  p_left  = NULL; }
  inline void                   Cluster::DeleteRight()    { delete p_right; p_right = NULL; }
  inline double                 Cluster::Mass2()    const { return m_momentum.Abs2(); }




  typedef std::list<Cluster*>          Cluster_List;
  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);


  /*!
    \file 
    \brief Contains the class AHADIC::Cluster and the corresponding struct AHADIC::ctp 
  */

  /*!
    \struct ctp
    \brief Short description of AHADIC::Cluster type. Available labels so far are: 
    none (0), quark-antiquark (11), quark-diquark (12), antidiquark-quark (21) antidiquark-diquark(22).
  */

  /*! 
    \class Cluster
    \brief The basic object of cluster fragmentation models: A colourless object composed of a
    quark-antiquark pair or similar (diquarks). In the realisation here, the objects are more or
    less carriers of information: flavour of the quarks, their momenta, links to their offsprings.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var ctp::code Cluster::m_type
    \brief The actual type of the cluster. In the default constructor it is set to "no", otherwise
    chosen corresponding to the flavours.
  */
  /*!
    \var ATOOLS::Flavour Cluster::m_flavours[2]
    \brief The flavours of the two constituents of the cluster. 
  */
  /*!
    \var ATOOLS::Vec4D Cluster::m_momentum
    \brief The momentum of the cluster.
  */
  /*!
    \var ATOOLS::Vec4D Cluster::m_momenta[2]
    \brief The momenta of the two constituents of the cluster. 
  */
  /*!
    \var ATOOLS::Poincare Cluster::m_boost
    \brief A boost operation connecting the lab frame with the cluster's rest frame.
  */
  /*!
    \var bool Cluster::m_hasboost
    \brief A flag indicating, whether a boost into the clusters rest rfame has already been initialised,
    set to "false" in the clusters constructors.
  */
  /*!
    \var Cluster * Cluster::p_left
    \brief A pointer to the left offspring (the offspring carrying the first flavour), set ot NULL in 
    the clusters constructors.
  */
  /*!
    \var Cluster * Cluster::p_right
    \brief A pointer to the right offspring (the offspring carrying the second flavour), set ot NULL in 
    the clusters constructors.
  */
  /*!
    \var Cluster * Cluster::p_prev
    \brief A pointer to the clusters predecessor, equals NULL if the cluster stems directly from the 
    cluster formation.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster::Cluster()
    \brief Default constructor, type is set to no, momenta to zero, and flavours to none.
  */
  /*!
    \fn Cluster::Cluster(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,const ATOOLS::Flavour &,const ATOOLS::Vec4D &)
    \brief Constructive constructor, the flavours are sorted such that the colour structure is triplet-antitriplet,
    and the type is then set accordingly. If no legal colour structure is found, the program aborts.
  */
  /*!
    \fn Cluster::~Cluster()
    \brief Destructor, deletes the pointers to the left and the right offspring.
  */
  /*!
    \fn Cluster * Cluster::GetLeft()  const
    \brief Returns the pointer to the left offspring cluster, p_left.
  */
  /*!
    \fn Cluster * Cluster::GetRight()  const
    \brief Returns the pointer to the right offspring cluster, p_right.
  */
  /*!
    \fn Cluster * Cluster::GetPrev()  const
    \brief Returns the pointer to the originator cluster, p_prev.
  */
  /*!
    \fn void Cluster::SetLeft(Cluster *);
    \brief Sets the pointer to the left offspring cluster, p_left. 
  */
  /*!
    \fn void Cluster::SetRight(Cluster *);
    \brief Sets the pointer to the right offspring cluster, p_right.
  */
  /*!
    \fn void Cluster::SetPrev(Cluster *);
    \brief Sets the pointer to the originator cluster, p_prev.
  */
  /*!
    \fn void Cluster::DeleteLeft();
    \brief Deletes the left offspring (and all its offsprings through its destructor) 
    and sets the pointer p_left to NULL.
  */
  /*!
    \fn void Cluster::DeleteRight();
    \brief Deletes the right offspring (and all its offsprings through its destructor) 
    and sets the pointer p_right to NULL.
  */
  /*!
    \fn ATOOLS::Flavour Cluster::GetFlav(const int) const
    \brief Returns the first or second flavour of m_flavours, the integer should be 1 or 2.
  */
  /*!
    \fn void Cluster::SetFlav(const int,const ATOOLS::Flavour &)
    \brief Sets the first or second flavour of m_flavours, the integer should be 1 or 2.
  */
  /*!
    \fn ATOOLS::Vec4D Cluster::Momentum(const int=0) const
    \brief Returns the momentum of one of the offsprings (argument equals 1 or 2) or the 
    clusters momentum (argument=0).
  */
  /*!
    \fn void Cluster::SetMomentum(const int,const ATOOLS::Vec4D &)
    \brief Sets the momentum of one of the offsprings (argument equals 1 or 2) or the 
    clusters momentum (argument=0).
  */
  /*!
    \fn double Cluster::Mass(const int=0) const;
    \brief Returns the constituent mass of one of the offsprings (argument equals 1 or 2) or 
    the clusters mass (argument=0), calculated from its momentum.
  */
  /*!
    \fn double Cluster::Mass2() const;
    \brief Returns the mass squared of the cluster, calculated from its momentum.
  */
  /*!
    \fn void Cluster::Update()
    \brief Refreshs the clusters momentum (equals then the sum of its two offsprings),
    checks for the colour structure of the constituents and sets the type accordingly. In
    case the colour configuration is illegal, the run aborts.
  */
  /*!
    \fn void Cluster::RescaleMomentum(ATOOLS::Vec4D)
    \brief Boosts the constituents into another system given by the the argument, which forms the
    clusters new momentum.

    In order to do this, first of all, the clusters actual c.m. system is constructed from its 
    actual momentum (which may be different from the sum of the two consituents). Then the
    constituents are boosted into this frame. From the argument a new c.m. system of the cluster 
    is constructed, the constituents are boosted back accordingly. Eventually, if the sum of
    the two new constituent momenta does not coincide with the new clusters momentum (the argument),
    a warning is printed to screen.
  */
  /*!
    \fn void Cluster::BoostInCMS()
    \brief Boosts the cluster, its constituents, and all offsprings into its c.m. system.

    If there is no boost so far (m_hasboost=false) the c.m. system of the cluster is constructed
    from its momentum. This momentum, the two constituents' momenta and all ofsprings are boosted into 
    the clusters' c.m. system. Then m_hasboost is set to true.    
  */
  /*!
    \fn void Cluster::BoostBack()
    \brief Boosts the cluster, its constituents, and all offsprings back into the lab system.

    If there is already a boost (m_hasboost=true) the clusters momentum, the two constituents' momenta 
    and all ofsprings are boosted back (into the lab system). Then m_hasboost is set to false.    
  */
  /*!
    \fn void Cluster::Boost(ATOOLS::Poincare &)
    \brief Boosts the clusters momentum, the momenta of its costituents and all offsprings by the
    given boost.
  */
  /*!
    \fn void Cluster::BoostBack(ATOOLS::Poincare &)
    \brief Boosts back the clusters momentum, the momenta of its costituents and all offsprings by the
    given boost.
  */
  /*!
    \fn void Cluster::BoostBack(ATOOLS::Vec4D &)
    \brief Uses the clusters boost to boost back a momentum.
  */
}
#endif
