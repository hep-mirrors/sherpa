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

  class Four_Fermion : public Cluster_Part {
  public:
    Four_Fermion() {}
    ~Four_Fermion() {}
    bool TestDecay(Cluster *,Part_List *) {}
  };
  /*!
    \file Cluster_Part.H
    \brief Contains the virtual mother class AHADIC::Cluster_Part and its explicit realisations
    AHADIC::QoverM and AHADIC::FourFermion. 
  */

  /*! 
    \class Cluster_Part
    \brief This class steers the cluster test decays into clusters, which may be kept or replaced by hadrons, if
    they are too light. This latter transition is managed by AHADIC::Hadron_Part.
  */
  /*! 
    \class QoverM
    \brief This is an explicit realisation of the virtual class AHADIC::Cluster_Part. Here, the kinematics of the
    test decays is fixed through simple momenta shuffling with a parameter Q.

    The kinematics of the cluster decay \f$C(p_1,p_2)\to C(q_1,q_2)+C(q_3,q_4)\f$ is given by
    \f[
    q_{1,4}^\mu = \left(1-\frac{Q}{M_C}\right)p_{1,2}^\mu \;\;\;\mbox{\rm and}\;\;\;
    q_{2,3}^\mu = \frac{Q}{M_C}p_{2,1}^\mu\,,
    \f]
    where \f$M_C^2=\sqrt{(p_1+p_2)^2}\f$ is the mass of the decaying cluster. The two new clusters thus have momentum
    \f[
    \begin{array}{lcl}
    Q_{12}^\mu &=& q_1^\mu+q_2^\mu = \left(1-\frac{Q}{M_C}\right)p_1^\mu + \frac{Q}{M_C}p_2^\mu\\
    Q_{34}^\mu &=& q_3^\mu+q_4^\mu = \left(1-\frac{Q}{M_C}\right)p_2^\mu + \frac{Q}{M_C}p_1^\mu
    \end{array}\,.
    \f]
  */
  /*! 
    \class Four_Fermion
    \brief This is an explicit realisation of the virtual class AHADIC::Cluster_Part. Here, the kinematics of the
    test decays is fixed through a matrix element for \f$S\to q\bar qq'\bar q'\f$. Not yet implemented.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Part
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Part::Cluster_Part()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Cluster_Part::~Cluster_Part()
    \brief Destructor, no task so far.
  */
  /*!
    \fn virtual bool Cluster_Part::TestDecay(Cluster *,Part_List *)
    \brief This method returns "true", if a test decay for the cluster in question according to the chosen explicit 
    strategy has successfully been constructed, otherwise "false" is returned.    
    
    Called by the Cluster_Decayer_Base::Treat(Cluster * cluster,Part_List * pl) method. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for QoverM
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Pair_Popper QoverM::m_popper
    \brief To pop out quark pairs in the test decays to find out whether the produced two objects fall into the
    cluster or the hadron regime.
  */
  /*!
    \var double QoverM::m_Q
    \brief The parameter for the cluster decays, at the moment \f$m_Q=1 \mbox{\rm GeV}^2\f$. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for QoverM
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn QoverM::QoverM()
    \brief Constructor, initialises the m_popper and m_Q.
  */
  /*!
    \fn QoverM::~QoverM()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool  QoverM::TestDecay(Cluster *,Part_List *)
    \brief This method produces a test decay of the cluster, including kinematics and a new quark (or diquark) pair;
    if such a test decay can be constructed the method yields "true", otherwise it yields "false".

    The cluster decay is performed in the decayers c.m. frame, kinematics of the cluster decay 
    \f$C(p_1,p_2)\to C(q_1,q_2)+C(q_3,q_4)\f$ is given by
    \f[
    q_{1,4}^\mu = \left(1-\frac{Q}{M_C}\right)p_{1,2}^\mu \;\;\;\mbox{\rm and}\;\;\;
    q_{2,3}^\mu = \frac{Q}{M_C}p_{2,1}^\mu\,,
    \f]
    where \f$M_C^2=\sqrt{(p_1+p_2)^2}\f$ is the mass of the decaying cluster. The two new clusters thus have momentum
    \f[
    \begin{array}{lcl}
    Q_{12}^\mu &=& q_1^\mu+q_2^\mu = \left(1-\frac{Q}{M_C}\right)p_1^\mu + \frac{Q}{M_C}p_2^\mu\\
    Q_{34}^\mu &=& q_3^\mu+q_4^\mu = \left(1-\frac{Q}{M_C}\right)p_2^\mu + \frac{Q}{M_C}p_1^\mu
    \end{array}\,.
    \f]
    Then a quark (or diquark) pair is determined with help of Pair_Popper::Pop(ATOOLS::Flavour &).
    A trivial check is performed whether the energies \f$E_{12}\f$ and \f$E_{34}\f$ of the new clusters
    are larger than the sum of the constituent masses. If not, a new trial flavour is popped; this is repeated until
    the conditions are met, then "true" is returned. If no flavour can be found that matches the fixed kinematics,
    "false" is returned.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Four_Fermion
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Four_Fermion::
    \brief 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Four_Fermion
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Four_Fermion::Four_Fermion()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Four_Fermion::~Four_Fermion()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool Four_Fermion::TestDecay(Cluster *,Part_List *)
    \brief 
  */
}
#endif
