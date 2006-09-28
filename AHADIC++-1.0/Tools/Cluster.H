#ifndef Cluster_H
#define Cluster_H

#include "Blob.H"
#include "Particle.H"
#include "Poincare.H"
#include "Message.H"
#include <list>
#include <iostream>

namespace AHADIC {
  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);

  struct control {
  public:
    static long int s_AHAparticles;
    static long int s_AHAprotoparticles;
    static long int s_AHAblobs;
  };

  struct ctp {
    enum code {
      no = 0,
      qq = 11,
      qd = 12,
      dq = 21,
      dd = 22
    };
  };

  struct ltp {
    enum code {
      no          = 0,
      leadingtrip = 1,
      leadinganti = 2,
      leadingboth = 3
    };
  };

  class Cluster {
  private:
    ctp::code           m_type;
    ltp::code           m_leads;
    ATOOLS::Flavour     m_flavours[2];
    ATOOLS::Vec4D       m_momentum, m_momenta[2];
    ATOOLS::Poincare    m_boost,m_rotate;
    bool                m_hasboost,m_hasrotate;

    Cluster           * p_left, * p_right, * p_prev;
    ATOOLS::Particle  * p_self;

    static long int     s_cluster_count, s_cluster_number;
    int                 m_number;
  public:
    Cluster();
    Cluster(const ATOOLS::Flavour &,const ATOOLS::Vec4D,
	    const ATOOLS::Flavour &,const ATOOLS::Vec4D);
    ~Cluster();

    ATOOLS::Blob     * CHHDecayBlob();

    Cluster          * GetLeft()  const;
    Cluster          * GetRight() const;
    Cluster          * GetPrev()  const;
    void               SetLeft(Cluster *);
    void               SetRight(Cluster *);
    void               SetPrev(Cluster *);
    ATOOLS::Particle * GetSelf() const;
    void               SetSelf(ATOOLS::Particle *);
    void               DeleteLeft();
    void               DeleteRight();
    ATOOLS::Flavour    GetFlav(const int) const;
    void               SetFlav(const int,const ATOOLS::Flavour &);
    ATOOLS::Vec4D      Momentum(const int=0) const;
    void               SetMomentum(const int,const ATOOLS::Vec4D);
    ltp::code          GetLeads() const;
    void               SetLeads(const ltp::code);
    double             Mass2() const;
    double             Mass(const int=0) const;
    void               Update();
    void               RescaleMomentum(ATOOLS::Vec4D);
    void               BoostInCMSAndRotateOnZ();
    void               RotateAndBoostBack();
    void               BoostInCMS();
    void               BoostBack();
    void               Boost(ATOOLS::Poincare &);
    void               BoostBack(ATOOLS::Poincare &);
    void               Rotate(ATOOLS::Poincare &);
    void               RotateBack(ATOOLS::Poincare &);
    void               BoostBack(ATOOLS::Vec4D &);
    void               RotateAndBoostBack(ATOOLS::Vec4D &);
    friend std::ostream& operator<<(std::ostream &,const Cluster &);

    inline int         Number()             { return m_number; }
    inline int         RemainingClusters()  { return s_cluster_count; }
    inline int         ClusterCounter()     { return s_cluster_number; }
    inline void        ResetClusterNumber() { s_cluster_number = 0; }
  };


  inline Cluster          * Cluster::GetLeft()                  const { return p_left;  }
  inline Cluster          * Cluster::GetRight()                 const { return p_right; }
  inline Cluster          * Cluster::GetPrev()                  const { return p_prev;  }
  inline void               Cluster::SetLeft(Cluster * left)          { p_left  = left;  }
  inline void               Cluster::SetRight(Cluster * right)        { p_right = right; }
  inline void               Cluster::SetPrev(Cluster * prev)          { p_prev  = prev;  }
  inline void               Cluster::DeleteLeft()  { delete p_left;  p_left  = NULL; }
  inline void               Cluster::DeleteRight() { delete p_right; p_right = NULL; }
  inline double             Cluster::Mass2()                    const { return m_momentum.Abs2(); }
  inline ltp::code          Cluster::GetLeads()                 const { return m_leads; }
  inline void               Cluster::SetLeads(const ltp::code l)      { m_leads = l; }
  inline ATOOLS::Particle * Cluster::GetSelf()                  const { return p_self; }
  inline void               Cluster::SetSelf(ATOOLS::Particle * part) { p_self = part; }




  typedef std::list<Cluster*>          Cluster_List;
  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);

  struct Proto_Particle {
  public:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_mom;
    char            m_info;
    Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info) :
      m_flav(flav), m_mom(mom), m_info(info) {}
  };

  typedef std::list<Proto_Particle * >        Proto_Particle_List;   
  typedef Proto_Particle_List::iterator       PPL_Iterator; 
  typedef Proto_Particle_List::const_iterator PPL_Const_Iterator; 
  typedef std::vector<Proto_Particle_List *>  VectorOfPPLs;
  typedef VectorOfPPLs::iterator              VPPL_Iterator;

  std::ostream & operator<<(std::ostream & s, const Proto_Particle_List & pl);
  /*!
    \file 
    \brief Contains the class AHADIC::Cluster and the structs AHADIC::ctp and
    AHADIC::ltp, which further descibe and label cluster properties, and the struct
    AHADIC::Proto_Particle.
  */

  /*!
    \struct ctp
    \brief Description of AHADIC::Cluster type. 

    Available labels are: 
    - none                = 0, 
    - quark-antiquark     = 11, 
    - quark-diquark       = 12, 
    - antidiquark-quark   = 21, 
    - antidiquark-diquark = 22.
    In practise, they are not used at the moment.
  */
  /*!
    \struct ltp
    \brief Identifying leading particles in a AHADIC::Cluster. 

    Available labels are: 
    - no          = 0,
    - leadingtrip = 1,
    - leadinganti = 2,
    - leadingboth = 3
    In practise, they are not used at the moment.
  */
  /*!
    \struct Proto_Particle 
    \brief This struct is used for the construction of clusters. 

    It acts as a container for the flavour and momentum of a particle and for
    the information, whether it is a leading particle or not.
  */
  /*! 
    \class Cluster
    \brief The basic object of cluster fragmentation models.

    The basic objects of this hadronisation model are colourless objects composed of 
    a quark-antiquark pair or similar (diquarks). In the realisation here, the objects 
    are more or less carriers of information: flavour of the quarks, their momenta, 
    links to their offsprings.

    During hadronisation, the clusters are formed (see AHADIC::Cluster_Formation_Handler
    and the classes called by it), after the gluons have been split into colour 
    triplet-antitriplet pairs, from quarks and antiquarks or corresponding diquark degrees
    of freedom.  Then the clusters are either transformed directly into hadrons
    (cluster->hadron and cluster->hadron+hadron transitions are handled in the
    AHADIC::Soft_Cluster_Handlder) or decayed into secondary clusters, which again may 
    be transformed into hadrons.  This stage of the hadronisation process is steered by
    the AHADIC::Cluster_Decay_Handler.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var ctp::code Cluster::m_type
    \brief The actual type of the cluster. 

    In the default constructor it is set to "no", otherwise chosen corresponding to the 
    flavours, cf. the description of AHADIC::ctp.
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
    \var ATOOLS::Poincare Cluster::m_rotate
    \brief A rotation, usually initialised such that the constituents of the cluster are 
    oriented along the z-axis.  
  */
  /*!
    \var bool Cluster::m_hasboost
    \brief A flag indicating, whether a boost into the clusters rest frame has 
    already been initialised, set to "false" in the clusters constructors.  This flag will
    be updated by each boost operation.  
  */
  /*!
    \var bool Cluster::m_hasrotate
    \brief A flag indicating, whether the rotation on the z-axis has already been initialised, 
    set to "false" in the clusters constructors.  This flag will be updated by each rotation. 
  */
  /*!
    \var Cluster * Cluster::p_left
    \brief A pointer to the left offspring (the offspring carrying the first flavour), set 
    to NULL in the clusters constructors.
  */
  /*!
    \var Cluster * Cluster::p_right
    \brief A pointer to the right offspring (the offspring carrying the second flavour), set 
    to NULL in the clusters constructors.
  */
  /*!
    \var Cluster * Cluster::p_prev
    \brief A pointer to the clusters predecessor, equals NULL if the cluster stems directly 
    from the cluster formation.
  */
  /*!
    \var ATOOLS::Particle  * p_self
    \brief A pointer to a particle representing the cluster - this can be used for extended
    event records.
  */
  /*!
    \var static long int s_cluster_count
    \brief For debugging reasons - a counter checking whether all clusters have been
    deleted.  
  */
  /*!
    \var static long int s_cluster_number
    \brief For debugging reasons - a quick way of identifying clusters.
  */
  /*!
    \var int m_number
    \brief For debugging reasons - a quick way of identifying clusters.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster::Cluster()
    \brief Default constructor 

    In this constructor, the type is set to 'no', all momenta are zero, and there
    are no flavours - hence ATOOLS::Flavour(ATOOLS::kf::none).  It is used only in the
    AHADIC::Soft_Cluster_Handler.  Also, the AHADIC::Cluster::s_cluster_count and
    AHADIC::Cluster::s_cluster_number are incremented and the cluster number is set to
    the actual (static) number.  
  */
  /*!
    \fn Cluster::Cluster(const ATOOLS::Flavour &,const ATOOLS::Vec4D,
                         const ATOOLS::Flavour &,const ATOOLS::Vec4D)
    \brief Constructive constructor

    Here, the flavours are sorted such that the colour structure is triplet-antitriplet,
    and the type is then set accordingly. If no legal colour structure is found, the 
    program aborts.  Also, momenta are set and the cluster's total momentum is set to
    the sum of its two constituents' momenta.  Also, the AHADIC::Cluster::s_cluster_count and
    AHADIC::Cluster::s_cluster_number are incremented and the cluster number is set to
    the actual (static) number.    
  */
  /*!
    \fn Cluster::~Cluster()
    \brief Destructor.

    In principle, the destructor doesn't do anything but decrementing the cluster counter
    AHADIC::Cluster::s_cluster_count.
  */
  /*!
    \fn ATOOLS::Blob * Cluster::CHHDecayBlob()
    \brief Constructs and returns a decay blob for the cluster.

    This method constructs a decay blob for the cluster with the type ATOOLS::btp::Cluster_Decay
    and the ATOOLS::blob_status::needs_hadrondecays.  The type specification, of course,
    is 'AHADIC-1.0'.  Then the cluster particle is added as incoming particle to the blob, and
    if existing, the offspring particles are added as outgoing particles.  If their flavour
    is not 'cluster', the corresponding offspring is deleted directly - it is a hadron now.
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
    \fn ATOOLS::Particle * Cluster::GetSelf() const
    \brief Returns the poniter to the cluster particle.
  */
  /*!
    \fn void Cluster::SetLeft(Cluster *);
    \brief Sets the pointer to the left offspring cluster, p_left. 

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::SetRight(Cluster *);
    \brief Sets the pointer to the right offspring cluster, p_right.

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::SetPrev(Cluster *);
    \brief Sets the pointer to the originator cluster, p_prev.

    This method is actually not being used in the code.
  */
  /*!
   \fn void Cluster::SetSelf(ATOOLS::Particle *)
   \brief Sets the cluster's particle pointer.
  */
  /*!
    \fn void Cluster::DeleteLeft();
    \brief Deletes the left offspring and sets the pointer p_left to NULL.

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::DeleteRight();
    \brief Deletes the right offspring and sets the pointer p_right to NULL.

    This method is actually not being used in the code.
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
    \fn void Cluster::SetMomentum(const int,const ATOOLS::Vec4D)
    \brief Sets the momentum of one of the offsprings (argument equals 1 or 2) or the 
    clusters momentum (argument=0).
  */
  /*!
    \fn ltp::code Cluster::GetLeads() const;
    \brief Returns the code for the leading particle configuration.
  */
  /*!
    \fn void Cluster::SetLeads(const ltp::code);
    \brief Sets the code for the leading particle configuration. 
  */
  /*!
    \fn double Cluster::Mass(const int=0) const;
    \brief Returns the constituent mass of one of the offsprings (argument equals 1 or 2) or 
    the clusters mass (argument=0), calculated directly from its momentum.
  */
  /*!
    \fn double Cluster::Mass2() const;
    \brief Returns the mass squared of the cluster, calculated directly from its momentum.
  */
  /*!
    \fn void Cluster::Update()
    \brief Performs on an actual cluster the same steps as in the constructive constructor, i.e.
    refreshs momenta and structure.

    This method refreshs the clusters momentum (equals then the sum of its two offsprings),
    checks for the colour structure of the constituents and sets the type accordingly. In
    case the colour configuration is illegal, the run aborts.  This method is used in
    AHADIC::Colour_Reconnections, where, eventually, constituents are exchanged between
    two clusters.
  */
  /*!
    \fn void Cluster::RescaleMomentum(ATOOLS::Vec4D)
    \brief Boosts the constituents into another system given by the the argument, which 
    forms the clusters new momentum.

    In order to do this, first of all, the clusters actual c.m. system is constructed from its 
    actual momentum (which may be different from the sum of the two consituents). Then the
    constituents are boosted into this frame. From the argument a new c.m. system of the cluster 
    is constructed, the constituents are boosted back accordingly. Eventually, if the sum of
    the two new constituent momenta does not coincide with the new clusters momentum (the argument),
    a warning is printed to screen.  This is the case, if the masses of the cluster related to
    the new and the old momentum are different.  
  */
  /*!
    \fn void Cluster::BoostInCMSAndRotateOnZ();
    \brief Boosts the cluster, its constituents, and all offsprings into its c.m. system and
    rotates it such that the constituents are parallel to the z-axis.  

    To this end, first Cluster::BoostInCMS() is called before a suitable rotation is constructed.
    It is defined such that the first constituent is oriented along the positive z-axis.  With
    this rotation Cluster::Rotate is called (and m_hasrotate is set to 'true').
  */
  /*!
    \fn void Cluster::RotateAndBoostBack();
    \brief Rotates and boosts back the cluster, its constituents, and all offsprings.

    For this the methods Cluster::RotateBack() and Cluster::BoostBack() are invoked. 
  */
  /*!
    \fn void Cluster::BoostInCMS()
    \brief Boosts the cluster, its constituents, and all offsprings into its c.m. system.

    If there is no boost so far (m_hasboost=false) the c.m. system of the cluster is 
    constructed from its momentum. This momentum, the two constituents' momenta and all 
    offsprings are boosted into the clusters' c.m. system. Finally, the cluster particle's 
    momentum is updated and m_hasboost is set to 'true'.    
  */
  /*!
    \fn void Cluster::BoostBack()
    \brief Boosts the cluster, its constituents, and all offsprings back into the lab system.

    If there is already a boost (m_hasboost=true) the clusters momentum, the two constituents' 
    momenta and all ofsprings are boosted back (into the lab system). Finally, the cluster particle's 
    momentum is updated and m_hasboost is set to 'false'.    
  */
  /*!
    \fn void Cluster::Boost(ATOOLS::Poincare &)
    \brief Boosts the clusters momentum, the momenta of its costituents and all offsprings by the
    given boost and updates the cluster particle's momentum correspondingly.
  */
  /*!
    \fn void Cluster::BoostBack(ATOOLS::Poincare &)
    \brief Boosts back the clusters momentum, the momenta of its costituents and all offsprings
    by the given boost and updates the cluster particle's momentum correspondingly.
  */
  /*!
    \fn void Cluster::Rotate(ATOOLS::Poincare &)
    \brief Rotates the clusters momentum, the momenta of its costituents and all offsprings by the
    given rotation and updates the cluster particle's momentum correspondingly.
  */
  /*!
    \fn void Cluster::RotateBack(ATOOLS::Poincare &)
    \brief Rotates back the clusters momentum, the momenta of its costituents and all offsprings 
    by the given rotation and updates the cluster particle's momentum correspondingly.
  */
  /*!
    \fn void Cluster::BoostBack(ATOOLS::Vec4D &)
    \brief Uses the clusters boost to boost back a momentum.
  */
  /*!
    \fn void Cluster::RotateAndBoostBack(ATOOLS::Vec4D &)
    \brief Uses the clusters boost and rotation to rotate and boost back a momentum.
  */
  /*!
    \fn friend std::ostream& operator<<(std::ostream &,const Cluster &)
    \brief Output the cluster in a stream.
  */

}
#endif
