#ifndef AHADIC_Tools_Cluster_H
#define AHADIC_Tools_Cluster_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include <list>
#include <iostream>

//#define AHAmomcheck
//#define memchecker

namespace AHADIC {
  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);

  class Proto_Particle;
  class Proto_Particle {
    static std::list<Proto_Particle * > s_actives;
  public:
    ATOOLS::Flavour  m_flav;
    ATOOLS::Vec4D    m_mom;
    char             m_info;
    double           m_mass, m_kt2max;
    Proto_Particle * p_partner;
  public:
    Proto_Particle(const Proto_Particle &);  
    Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info);  
    ~Proto_Particle();

    bool CheckConsistency(std::ostream & s,std::string method);

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Proto_Particles -----"<<std::endl;
      for (std::list<Proto_Particle * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->m_flav<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };
  std::ostream & operator<<(std::ostream & s, const Proto_Particle &);



  class Proto_Particle_List;
  class Proto_Particle_List : public std::list<Proto_Particle * > {
  private:
    static std::list<Proto_Particle_List * > s_actives;
  public:
    Proto_Particle_List() {
      s_actives.push_back(this);
    }

    Proto_Particle_List(const Proto_Particle_List & pl) {
      s_actives.push_back(this);
      for (std::list<Proto_Particle * >::const_iterator piter=pl.begin();
	   piter!=pl.end();piter++)
	push_back((*piter));
    }

    ~Proto_Particle_List() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Proto_Particle_Lists -----"<<std::endl;
      for (std::list<Proto_Particle_List * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef Proto_Particle_List::iterator       PPL_Iterator; 
  typedef Proto_Particle_List::const_iterator PPL_Const_Iterator; 
  std::ostream & operator<<(std::ostream & s, const Proto_Particle_List & pl);


  class ListOfPPLs;
  class ListOfPPLs : public std::list<Proto_Particle_List * > {
  private:
    static std::list<ListOfPPLs * > s_actives;
  public:
    ListOfPPLs() {
      s_actives.push_back(this);
    }

    ~ListOfPPLs() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active ListOfPPLss -----"<<std::endl;
      for (std::list<ListOfPPLs * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef ListOfPPLs::iterator                LPPL_Iterator;

  struct control {
  public:
    static long int s_AHAparticles;
    static long int s_AHAprotoparticles;
    static long int s_AHAblobs;
  };

  class Cluster;

  class Cluster {
  private:
    bool                m_active;
    Proto_Particle    * p_trip, * p_anti;
    ATOOLS::Vec4D       m_momentum;
    ATOOLS::Flavour     m_flav;
    double              m_kt, m_z;

    ATOOLS::Poincare    m_boost,m_rotate;
    bool                m_hasboost,m_hasrotate;

    Cluster           * p_left, * p_right;
    Cluster           * p_prev;
    Cluster           * p_nbtrip, * p_nbanti;
    long int            m_number;

    std::vector<ATOOLS::Flavour> m_decayproducts;

    static long int     s_cluster_count, s_cluster_number;
    static std::list<Cluster *> s_actives;

    void               EnsureMomentum();
  public:
    Cluster(ATOOLS::Vec4D=ATOOLS::Vec4D(0.,0.,0.,0.),
	    ATOOLS::Flavour=ATOOLS::Flavour(kf_cluster),
	    bool active=true);
    Cluster(Proto_Particle *,Proto_Particle *);
    ~Cluster();

    void               Update();
    bool               CheckConsistency(std::ostream & s,std::string method);
    ATOOLS::Blob     * ConstructDecayBlob();

    ATOOLS::Particle * GetSelf() const;

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

    bool               Active() const;
    ATOOLS::Flavour    GetFlav() const;
    Cluster *          GetLeft()  const;
    Cluster *          GetRight() const;
    Cluster *          GetPrev()  const;
    Cluster *          GetNBTrip()  const;
    Cluster *          GetNBAnti() const;
    void               SetActive(bool);
    void               SetFlav(ATOOLS::Flavour);
    void               SetLeft(Cluster *);
    void               SetRight(Cluster *);
    void               SetPrev(Cluster *);
    void               SetNBTrip(Cluster *);
    void               SetNBAnti(Cluster *);


    Proto_Particle *   GetTrip() const;
    Proto_Particle *   GetAnti() const;
    void               SetTrip(Proto_Particle *);
    void               SetAnti(Proto_Particle *);
    ATOOLS::Vec4D      Momentum() const;
    void               SetMomentum(const ATOOLS::Vec4D);
    double             Mass2() const;
    double             Mass() const;
    double             KT() const;
    double             Z() const;
    void               SetKT(const double);
    void               SetZ(const double);

    void               push_back(ATOOLS::Flavour);
    void               clear();
    size_t             size() const;
    ATOOLS::Flavour    operator[](size_t pos);


    void               Print();
    void               Delete();

    friend std::ostream& operator<<(std::ostream &,const Cluster &);

    inline long int    Number()             { return m_number; }
    static long int    RemainingClusters()  { return s_cluster_count; }
    static void        ResetClusterNumber() { s_cluster_number = 0; }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Clusters -----"<<std::endl;
      for (std::list<Cluster * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->m_number<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };


  
  inline bool            Cluster::Active() const              { return m_active;  }
  inline ATOOLS::Flavour Cluster::GetFlav() const             { return m_flav;    }
  inline Cluster *       Cluster::GetLeft() const             { return p_left;    }
  inline Cluster *       Cluster::GetRight() const            { return p_right;   }
  inline Cluster *       Cluster::GetPrev() const             { return p_prev;    }
  inline Cluster *       Cluster::GetNBTrip() const           { return p_nbtrip;  }
  inline Cluster *       Cluster::GetNBAnti() const           { return p_nbanti;  }
  inline void            Cluster::SetActive(bool act)         { m_active = act;   }
  inline void            Cluster::SetFlav(ATOOLS::Flavour fl) { m_flav   = fl;    }
  inline void            Cluster::SetLeft(Cluster * left)     { p_left   = left;  }
  inline void            Cluster::SetRight(Cluster * right)   { p_right  = right; }
  inline void            Cluster::SetPrev(Cluster * prev)     { p_prev   = &*prev;}
  inline void            Cluster::SetNBTrip(Cluster * clu)    { p_nbtrip = clu;   }
  inline void            Cluster::SetNBAnti(Cluster * clu)    { p_nbanti = clu;   }


  inline Proto_Particle * Cluster::GetTrip() const                { return p_trip; }
  inline Proto_Particle * Cluster::GetAnti() const                { return p_anti; }
  inline void             Cluster::SetTrip(Proto_Particle * trip) { p_trip = trip; }
  inline void             Cluster::SetAnti(Proto_Particle * anti) { p_anti = anti; }
  inline ATOOLS::Vec4D    Cluster::Momentum() const                   { return m_momentum; }
  inline void             Cluster::SetMomentum(const ATOOLS::Vec4D p) { m_momentum = p; }
  inline double           Cluster::Mass2() const             { return m_momentum.Abs2(); }
  inline double           Cluster::Mass() const              { return sqrt(m_momentum.Abs2()); }
  inline void             Cluster::push_back(ATOOLS::Flavour flav) { m_decayproducts.push_back(flav); }
  inline void             Cluster::clear()                         { m_decayproducts.clear(); }
  inline size_t           Cluster::size() const                    { return m_decayproducts.size(); }
  inline ATOOLS::Flavour  Cluster::operator[](size_t pos)          { return m_decayproducts[pos]; }
  inline double           Cluster::KT() const                      { return m_kt; }
  inline double           Cluster::Z() const                       { return m_z; }
  inline void             Cluster::SetKT(const double kt)          { m_kt = kt; }
  inline void             Cluster::SetZ(const double z)            { m_z = z; }


  class Cluster_List;
  class Cluster_List : public std::list<Cluster * > {
  private:
    static std::list<Cluster_List * > s_actives;
  public:
    Cluster_List() {
      s_actives.push_back(this);
    }

    ~Cluster_List() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Cluster_Lists -----"<<std::endl;
      for (std::list<Cluster_List * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);
  /*!
    \file 
    \brief Contains the class AHADIC::Cluster, the structs AHADIC::control and 
    AHADIC::Proto_Particle.
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
    are more or less carriers of information: flavour of the quarks and their momenta,
    both encoded in a Proto_Particle, links to the clusters' offsprings.

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
  // Description of member variables for struct Proto_Particle
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var ATOOLS::Flavour Proto_Particle::m_flav
    \brief The original particles flavour, copied in this little container.
  */
  /*!
    \var ATOOLS::Vec4D Proto_Particle::m_mom
    \brief The original particles four-momentum, copied in this little container.
  */
  /*!
    \var char Proto_Particle::m_info
    \brief The original particles info-char, copied in this little container.
  */
  /*!
    \fn Proto_Particle::Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info) 
    \brief Explicit constructor, increments the counter::AHAprotoparticles counter.
  */
  /*!
    \fn Proto_Particle::Proto_Particle(const Proto_Particle &) 
    \brief Copy constructor, increments the counter::AHAprotoparticles counter.
  */
  /*!
    \fn Proto_Particle::~Proto_Particle()
    \brief Destructor, decrements the counter::AHAprotoparticles counter.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Proto_Particle * Cluster::p_trip
    \brief The Proto_Particle with triplet colour charge (quark or antidiquark).
  */
  /*!
    \var Proto_Particle * Cluster::p_antitrip
    \brief The Proto_Particle with antitriplet colour charge (antiquark or diquark).
  */
  /*!
    \var ATOOLS::Vec4D Cluster::m_momentum
    \brief The momentum of the cluster.
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
    \var ATOOLS::Particle * Cluster::p_self
    \brief A pointer to a particle representing the cluster - this can be used for extended
    event records.
  */

  /*!
    \var static long int Cluster::s_cluster_count
    \brief For debugging reasons - a counter checking whether all clusters have been
    deleted.  
  */
  /*!
    \var static long int Cluster::s_cluster_number
    \brief For debugging reasons - a quick way of identifying clusters.
  */
  /*!
    \var long int Cluster::m_number
    \brief For debugging reasons - a quick way of identifying clusters.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster::Cluster()
    \brief Default constructor 

    In this constructor, the type is set to 'no', all momenta are zero, and there
    are no flavours - hence ATOOLS::Flavour(kf_none).  It is used only in the
    AHADIC::Soft_Cluster_Handler.  Also, the AHADIC::Cluster::s_cluster_count and
    AHADIC::Cluster::s_cluster_number are incremented and the cluster number is set to
    the actual (static) number.  
  */
  /*!
    \fn Cluster::Cluster(Proto_Particle * trip, Proto_Particle * anti)
    \brief Constructive constructor

    Here, the Proto_Particle's must be sorted such that the colour structure is 
    triplet-antitriplet.  If no legal colour structure is found, the program aborts.  
    The cluster's total momentum is set to the sum of its two constituents' momenta.  Also, 
    the AHADIC::Cluster::s_cluster_count and AHADIC::Cluster::s_cluster_number are incremented 
    and the cluster number is set to the actual (static) number.    
  */
  /*!
    \fn Cluster::~Cluster()
    \brief Destructor.

    The destructor deletes the two Proto_Particle's and decrements the cluster counter 
    AHADIC::Cluster::s_cluster_count.
  */
  /*!
    \fn ATOOLS::Blob * Cluster::DecayBlob()
    \brief Constructs and returns a decay blob for the cluster.

    This method constructs a decay blob for the cluster with the type ATOOLS::btp::Cluster_Decay
    and the ATOOLS::blob_status::needs_hadrondecays.  The type specification, of course,
    is 'AHADIC-1.0'.  Then the cluster particle is added as incoming particle to the blob, and
    if existing, the offspring particles are added as outgoing particles.  If their flavour
    is not 'cluster', the corresponding offspring is deleted directly - it is a hadron now.
  */
  /*!
    \fn Cluster * Cluster::GetLeft()  const
    \brief Returns the pointer to the left offspring cluster, Cluster::p_left.
  */
  /*!
    \fn Cluster * Cluster::GetRight()  const
    \brief Returns the pointer to the right offspring cluster, Cluster::p_right.
  */
  /*!
    \fn Cluster * Cluster::GetPrev()  const
    \brief Returns the pointer to the originator cluster, Cluster::p_prev.
  */
  /*!
    \fn void Cluster::SetLeft(Cluster *);
    \brief Sets the pointer to the left offspring cluster, Cluster::p_left. 

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::SetRight(Cluster *);
    \brief Sets the pointer to the right offspring cluster, Cluster::p_right.

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::SetPrev(Cluster *);
    \brief Sets the pointer to the originator cluster, Cluster::p_prev.

    This method is actually not being used in the code.
  */
  /*!
    \fn ATOOLS::Particle * Cluster::GetSelf() const
    \brief Returns the poniter to the cluster particle.
  */
  /*!
   \fn void Cluster::SetSelf(ATOOLS::Particle *)
   \brief Sets the cluster's particle pointer.
  */
  /*!
    \fn Proto_Particle * Cluster::GetTrip()  const
    \brief Returns the pointer to the triplet constituent, Cluster::p_trip.
  */
  /*!
    \fn Proto_Particle * Cluster::GetAnti()  const
    \brief Returns the pointer to the antitriplet constituent, Cluster::p_anti.
  */
  /*!
    \fn void Cluster::SetTrip(Proto_Particle *);
    \brief Sets the pointer to the triplet constituent, Cluster::p_trip. 
  */
  /*!
    \fn void Cluster::SetAnti(Proto_Particle *);
    \brief Sets the pointer to the antitriplet constituent, Cluster::p_anti.
  */

  /*!
    \fn void Cluster::DeleteLeft();
    \brief Deletes the left offspring and sets the pointer Cluster::p_left to NULL.

    This method is actually not being used in the code.
  */
  /*!
    \fn void Cluster::DeleteRight();
    \brief Deletes the right offspring and sets the pointer Cluster::p_right to NULL.

    This method is actually not being used in the code.
  */
  /*!
    \fn ATOOLS::Vec4D Cluster::Momentum() const
    \brief Returns the clusters momentum.
  */
  /*!
    \fn void Cluster::SetMomentum()
    \brief Sets the clusters momentum.
  */
  /*!
    \fn double Cluster::Mass() const;
    \brief Returns the clusters mass, calculated directly from its momentum.
  */
  /*!
    \fn double Cluster::Mass2() const;
    \brief Returns the mass squared of the cluster, calculated directly from its momentum.
  */
  /*!
    \fn void Cluster::Update()
    \brief Performs on an actual cluster the same steps as in the constructive constructor, 
    i.e. refreshs the momentum and checks the structure.
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
    this rotation Cluster::Rotate is called (and Cluster::m_hasrotate is set to 'true').
  */
  /*!
    \fn void Cluster::RotateAndBoostBack();
    \brief Rotates and boosts back the cluster, its constituents, and all offsprings.

    For this the methods Cluster::RotateBack() and Cluster::BoostBack() are invoked. 
  */
  /*!
    \fn void Cluster::BoostInCMS()
    \brief Boosts the cluster, its constituents, and all offsprings into its c.m. system.

    If there is no boost so far (Cluster::m_hasboost=false) the c.m. system of the cluster is 
    constructed from its momentum. This momentum, the two constituents' momenta and all 
    offsprings are boosted into the clusters' c.m. system. Finally, the cluster particle's 
    momentum is updated and Cluster::m_hasboost is set to 'true'.    
  */
  /*!
    \fn void Cluster::BoostBack()
    \brief Boosts the cluster, its constituents, and all offsprings back into the lab system.

    If there is already a boost (Cluster::m_hasboost=true) the clusters momentum, the two 
    constituents' momenta and all ofsprings are boosted back (into the lab system). Finally, 
    the cluster particle's momentum is updated and Cluster::m_hasboost is set to 'false'.    
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
  /*!
    \fn inline int Cluster::Number()
    Returns the identifier number of the cluster.  For debugging only.
  */
  /*!
    \fn inline int RemainingClusters()
    Returns the number of undeleted clusters.  For debugging only.
  */
  /*!
    \fn inline void ResetClusterNumber()
    Resets the counter.  For debugging only.
  */
}
#endif
