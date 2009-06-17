#ifndef AHADIC_Decays_Cluster_Part_H
#define AHADIC_Decays_Cluster_Part_H


#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Dipole_Splitter.H"
#include "ATOOLS/Math/Histogram.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Part {
  protected:
    bool              m_ana;
    double            m_pt2max_factor;
    Dipole_Splitter * p_splitter;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    double PT2Max(Cluster * cluster) const;
  public:
    Cluster_Part(Dipole_Splitter * =NULL,bool=false);
    ~Cluster_Part();
    bool TestDecay(Cluster * const);
  };

  /*!
    \file Cluster_Part.H
    \brief This file contains the class Cluster_Part.
  */
  /*!
    \class Cluster_Part
    \brief The class cluster part is responsible for the decays of clusters into test clusters.
    These are either kept as clusters for further decays or they are transformed into hadrons by
    the Soft_Cluster_Handler.
  */
  /*!
    \var int Cluster_Part::m_4Qmode
    \brief Up to now: 0 by default. 

    \todo: Think about improvement.
  */
  /*!
    \var int Cluster_Part::m_ybar_mode
    \brief  Up to now: 0 by default.  For details, cf. 
    Cluster_Part::GetYBar(const double,const double)
  */
  /*!
    \var double Cluster_Part::m_fraction
    \brief The variable \f$f\f$ in Cluster_Part::BuildKinematics(Cluster * const,ATOOLS::Flavour &).

    It is set by the parameter "CLUSTER_MASS_FRACTION", default is 0.5.
  */
  /*!
    \var double Cluster_Part::m_ystar_expvalue
    \brief The expectation value of the \f$y^*\f$ shift in terms of \f$y^{*,max}\f$.

    For more details see the description of 
    Cluster_Part::GetYStar(const double,const double), there it is named 
    \f$\langle\delta y^*\rangle\f$.
    It is set by the parameter "<Y*_SHIFT>", default is 0.5.
  */
  /*!
    \var double Cluster_Part::m_ystar_sigma
    \brief The width of the \f$y^*\f$ distribution in a Gaussian.

    For more details see the description of 
    Cluster_Part::GetYStar(const double,const double), there it is named  \f$\sigma_y\f$.
    It is set by the parameter "<Y*_WIDTH>", default is 0.5.
  */
  /*!
    \var ATOOLS::Vec4D Cluster_Part::m_momenta[4]
    \brief The four momenta of the cluster constituents and the popped quark pair.
  */
  /*!
    \var bool Cluster_Part::m_ana
    \brief Switches the analysis on and off.  Default is true.
   */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> Cluster_Part::m_histograms
    \brief A number of histograms monitoring the decay kinematics.

    In particular, the \f$p_\perp\f$ distribution and the centre-of-mass energy squared of
    the popped quarks is monitored, as well as the flavour composition of the popped quarks.
    In addition, \f$y^*\f$, \f$y^{*}/y^{*, max}\f$, \f$\bar y\f$, \f$\bar y/\bar y^{max}\f$
    are monitored.  
  */
  /*!
    \fn Cluster_Part::Cluster_Part()
    \brief Initialises the constants and, if neccessary, the histograms.
  */
  /*!
    \fn Cluster_Part::~Cluster_Part()
    \brief Outputs the histograms from the analysis.
  */
  /*!
    \fn bool Cluster_Part::TestDecay(Cluster *)
    \brief Decides the popped flavour and instantiates two clusters accordingly.

    The procedure is as follows:
    - First the decaying cluster is boosted in its c.m. frame and its constitutents
      are oriented along the z-axis;
    - then Cluster_Part::BuildKinematics(Cluster * const,ATOOLS::Flavour &) is invoked
      to pop a flavour and fix the decay kinematics;
    - finally, offspring clusters are constructed with flavour content and momenta
      (eventually these clusters may become hadrons through the Soft_Cluster_Handler);
    - the system is boosted and rotated back, the offsprings are linked to the decaying 
      cluster, corresponding particles are instantiated for them.
  */
  /*!
    \fn bool Cluster_Part::UpdateDecay(Cluster *,const int)
    \brief Updates decay kinematics after one or both outgoing particles have been
    shifted to a new mass-shell.

    The basic idea is that the original decay is boosted in its c.m. system, with orientation
    of the original clusters constituents along the z-axis.  The energies and longitudinal (z-) 
    momenta of the outgoing particles are then recalculated with the new shifted masses of
    the offsrpings.  Their \f$p_\perp\f$ is kept and eventually rotated around the z-axis.
    Finally the full system is boosted and rotated back in its original frame.
  */
  /*!
    \fn bool Cluster_Part::BuildKinematics(Cluster * const,ATOOLS::Flavour &)
    \brief Builds the kinematics of a cluster->cluster+cluster decay in its own
    rest frame.

    The algorithm is as follows:
      The procedure is repeated, until physical energies for the original cluster
      constituents (larger than their mass) emerge.  

  */
}
#endif
