#ifndef Cluster_Part_H
#define Cluster_Part_H

#include "Pair_Popper.H"
#include "Cluster.H"
#include "Histogram.H"
#include <map>
#include "Return_Value.H"


namespace AHADIC {
  class Cluster_Part {
  protected:
    Pair_Popper m_popper;
    int         m_4Qmode;
    bool        m_ana;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;
    bool Veto(Cluster *,const ATOOLS::Flavour &);

    virtual bool BuildKinematics(Cluster *,ATOOLS::Flavour &,
				 ATOOLS::Vec4D *,double *)=0;
  public:
    Cluster_Part(const double,const int);
    virtual ~Cluster_Part() {};
    ATOOLS::Return_Value::code TestDecay(Cluster *);

    inline virtual void SetDecayAnalysisOn() { m_ana=true; }
  };

  class Schwinger_Uniform : public Cluster_Part {
    double m_tension, m_fraction, m_ystar_expvalue,m_ystar_sigma;

    double GetYStar(const double,const double);
    bool   BuildKinematics(Cluster *,ATOOLS::Flavour &,
			   ATOOLS::Vec4D *,double *);
  public:
    Schwinger_Uniform();
    ~Schwinger_Uniform();

    void SetDecayAnalysisOn();
  };

  class Four_Fermion : public Cluster_Part {
    bool BuildKinematics(Cluster *,ATOOLS::Flavour &,
			 ATOOLS::Vec4D *,double *) { abort(); }
  public:
    Four_Fermion() : Cluster_Part(0.,0) {}
    ~Four_Fermion() {}
  };

  class Cluster_Part_QBreakup : public Cluster_Part {
  protected:
    double m_Q;
    virtual const double SelectQ(const double) const=0;
    bool BuildKinematics(Cluster *,ATOOLS::Flavour &,
			 ATOOLS::Vec4D *,double *); 
  public:
    inline Cluster_Part_QBreakup(const double tension=0.,const int mode=1) : 
      Cluster_Part(tension,mode),
      m_Q(hadpars.Get(std::string("Q_breakup"))) {}
    virtual ~Cluster_Part_QBreakup() {}
  };

  class Simple_Q_over_M : public Cluster_Part_QBreakup {
    const double SelectQ(const double) const { return m_Q; }
  public:
    Simple_Q_over_M();
    ~Simple_Q_over_M() {}
  };

  class Running_Q_over_M : public Cluster_Part_QBreakup {
    const double SelectQ(const double) const;
  public:
    Running_Q_over_M();
    ~Running_Q_over_M() {}
  };

  /*!
    \file Cluster_Part.H
    \brief Contains the virtual mother class AHADIC::Cluster_Part and its explicit realisations
    AHADIC::Simple_Q_over_M, AHADIC::Running_Q_over_M, and AHADIC::FourFermion. 
  */

  /*! 
    \class Cluster_Part
    \brief This class steers the cluster test decays into clusters, which may be kept or replaced by hadrons, if
    they are too light. This latter transition is managed by AHADIC::Hadron_Part.
  */
  /*! 
    \class Simple_Q_over_M
    \brief This is an explicit realisation of the virtual class AHADIC::Cluster_Part. Here, the kinematics of the
    test decays is fixed through simple momenta shuffling with a fixed parameter Q.

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
    \class Running_Q_over_M
    \brief This is an explicit realisation of the virtual class AHADIC::Cluster_Part. Here, the kinematics of the
    test decays is fixed through simple momenta shuffling with a parameter varying Q.

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
    The parameter Q is selected according to
    \f[
    {\cal P}(Q) \sim \exp\left[-\frac{Q^2}{M_C^2}\right]
    \f]
    in the interval
    \f[
    Q\in\left[Q_0,\,Q_0+\sqrt{Q_0M_C}\right]
    \f]
    with a free parameter \f$Q_0\f$.
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
  // Description of member variables for Simple_Q_over_M
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Pair_Popper Simple_Q_over_M::m_popper
    \brief To pop out quark pairs in the test decays to find out whether the produced two objects fall into the
    cluster or the hadron regime.
  */
  /*!
    \var double Simple_Q_over_M::m_Q
    \brief The parameter for the cluster decays, at the moment \f$m_Q=1 \mbox{\rm GeV}^2\f$. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Simple_Q_over_M
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Simple_Q_over_M::Simple_Q_over_M()
    \brief Constructor, initialises the m_popper and m_Q.
  */
  /*!
    \fn Simple_Q_over_M::~Simple_Q_over_M()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool  Simple_Q_over_M::TestDecay(Cluster *,Part_List *)
    \brief This method produces a test decay of the cluster, including kinematics and a new quark (or diquark) pair;
    if such a test decay can be constructed the method yields "true", otherwise it yields "false".

    The cluster decay is performed in the decayers c.m. frame, kinematics of the cluster decay 
    \f$C(p_1,p_2)\to C(q_1,q_2)+C(q_3,q_4)\f$ is given by
    \f[
    q_{1,4}^\mu = \left(1-\frac{Q_0}{M_C}\right)p_{1,2}^\mu \;\;\;\mbox{\rm and}\;\;\;
    q_{2,3}^\mu = \frac{Q_0}{M_C}p_{2,1}^\mu\,,
    \f]
    where \f$M_C^2=\sqrt{(p_1+p_2)^2}\f$ is the mass of the decaying cluster. The two new clusters thus have momentum
    \f[
    \begin{array}{lcl}
    Q_{12}^\mu &=& q_1^\mu+q_2^\mu = \left(1-\frac{Q_0}{M_C}\right)p_1^\mu + \frac{Q_0}{M_C}p_2^\mu\\
    Q_{34}^\mu &=& q_3^\mu+q_4^\mu = \left(1-\frac{Q_0}{M_C}\right)p_2^\mu + \frac{Q_0}{M_C}p_1^\mu
    \end{array}\,.
    \f]
    This leads to masses of, roughly,
    \f[
    M_{12}^2 = Q_{12}^2 = m_1^2\left(1-\frac{2Q_0}{M_C}\right) + \frac{Q_0^2}{M_C^2}(m_1^2+m_2^2) - 
                          2(p_1\cdot p_2)\left(\frac{Q_0}{M_C}-\frac{Q_0^2}{M_C^2}\right) \approx Q_0M_C\,,
    \f]
    i.e. a significant reduction for the emerging clusters w.r.t. the original one. <br>
    Then a quark (or diquark) pair is determined with help of Pair_Popper::Pop(ATOOLS::Flavour &).
    A trivial check is performed whether the energies \f$E_{12}\f$ and \f$E_{34}\f$ of the new clusters
    are larger than the sum of the constituent masses. If not, a new trial flavour is popped; this is repeated until
    the conditions are met, then "true" is returned. If no flavour can be found that matches the fixed kinematics,
    "false" is returned.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Running_Q_over_M
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Pair_Popper Running_Q_over_M::m_popper
    \brief To pop out quark pairs in the test decays to find out whether the produced two objects fall into the
    cluster or the hadron regime.
  */
  /*!
    \var double Running_Q_over_M::m_Q
    \brief The parameter for the cluster decays, at the moment \f$m_Q=1 \mbox{\rm GeV}^2\f$. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Running_Q_over_M
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Running_Q_over_M::Running_Q_over_M()
    \brief Constructor, initialises the m_popper and m_Q.
  */
  /*!
    \fn Running_Q_over_M::~Running_Q_over_M()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool  Running_Q_over_M::TestDecay(Cluster *,Part_List *)
    \brief This method produces a test decay of the cluster, including kinematics and a new quark (or diquark) pair;
    if such a test decay can be constructed the method yields "true", otherwise it yields "false".

    Having selected an actual value \f$Q\f$ with help of SelectQ(const double), the cluster decay is performed in the 
    decayers c.m. frame, kinematics of the cluster decay 
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
    This leads to masses of, roughly,
    \f[
    M_{12}^2 = Q_{12}^2 = m_1^2\left(1-\frac{2Q}{M_C}\right) + \frac{Q^2}{M_C^2}(m_1^2+m_2^2) - 
                          2(p_1\cdot p_2)\left(\frac{Q}{M_C}-\frac{Q^2}{M_C^2}\right) \approx QM_C\,,
    \f]
    i.e. a significant reduction for the emerging clusters w.r.t. the original one. <br>
    Then a quark (or diquark) pair is determined with help of Pair_Popper::Pop(ATOOLS::Flavour &).
    A trivial check is performed whether the energies \f$E_{12}\f$ and \f$E_{34}\f$ of the new clusters
    are larger than the sum of the constituent masses. If not, a new trial flavour is popped; this is repeated until
    the conditions are met, then "true" is returned. If no flavour can be found that matches the fixed kinematics,
    "false" is returned.
  */
  /*!
    \fn double Running_Q_over_M::SelectQ(const double)
    \brief Given the cluster's mass (the argument), an actual value \f$Q\f$ is selected according to the 
    probability density
    \f[
    {\cal P}(Q) \sim \exp\left[-\frac{Q^2}{M_C^2}\right]
    \f]
    in the interval
    \f[
    Q\in\left[Q_0,\,Q_0+\sqrt{Q_0M_C}\right]
    \f]
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
