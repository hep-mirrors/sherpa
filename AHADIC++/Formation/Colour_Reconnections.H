#ifndef AHADIC_Formation_Colour_Reconnections_H
#define AHADIC_Formation_Colour_Reconnections_H

#include "AHADIC++/Tools/Cluster.H"

namespace AHADIC {
  class Colour_Reconnections {
  private:
    int    m_kinmode, m_colmode;
    double m_t0,m_w12,m_w34,m_w14,m_w23;

    bool   TestClusters(Cluster * cl1,Cluster * cl2,int gen=1);
    double KinematicWeight(const ATOOLS::Vec4D & mom1,const ATOOLS::Vec4D & mom2,
			   const ATOOLS::Vec4D & mom3,const ATOOLS::Vec4D & mom4);
    double ColourWeight(int gen=1);
  public:
    Colour_Reconnections(int kinmode=1,int colmode=1,double t0=1.);
    ~Colour_Reconnections();
    void   Singlet_CR(Cluster_List * clin);
    void   Two_Singlet_CR(Cluster_List * cl1,Cluster_List * cl2);
  };
  /*!
    \file Colour_Reconnections.H
    \brief Contains the class AHADIC::Colour_Reconnections.
  */
  /*! 
    \class Colour_Reconnections
    \brief Models non-pertubative colour reconnections that lead to a reshuffling of the colour
    flow in individual events.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Colour_Reconnections
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var int Colour_Reconnections::m_kinmode
    \brief The mode for the selection of the kinematical weight.

    Cases:
    - m_kinmode = 0: Kinematical weight is always zero, hence no reshuffling;
    - m_kinmode = 1: Kinematical weight is constructed from invariant masses (default value);
    - m_kinmode = 2: Kinematical weight is constructed from relative transverse momenta (foreseen).
  */
  /*!
    \var int Colour_Reconnections::m_colmode
    \brief The mode for the selection of the colour weight.

    Cases:
    - m_colmode = 0: Colour weight equals zero,  hence no reshuffling;
    - m_colmode = 1: Colour weight equals 1/9.
  */
  /*!
    \var double Colour_Reconnections::m_t0
    \brief The offset fro the construction of the kinematical weight, default is \f$1 \mbox{\rm GeV}^2\f$.
  */
  /*!
    \var double Colour_Reconnections::m_w12
    \brief The kinematical weight \f$w_{12}\f$, see below,
    KinematicWeight(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &).
 */
  /*!
    \var double Colour_Reconnections::m_w34
    \brief The kinematical weight \f$w_{34}\f$, see below,
    KinematicWeight(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &).
  */
  /*!
    \var double Colour_Reconnections::m_w14
    \brief The kinematical weight \f$w_{14}\f$, see below,
    KinematicWeight(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &).
  */
  /*!
    \var double Colour_Reconnections::m_w23
    \brief The kinematical weight \f$w_{23}\f$, see below,
    KinematicWeight(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Colour_Reconnections
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Colour_Reconnections::Colour_Reconnections(int kinmode=1,int colmode=1,double t0=1.)
    \brief Constructor, sets the parameters for the colour reconnections.
  */
  /*!
    \fn Colour_Reconnections::~Colour_Reconnections()
    \brief Destructor, so far it has no task.
  */
  /*!
    \fn void Colour_Reconnections::Singlet_CR(Cluster_List * clin)
    \brief Models colour reconnections inside one singlet list.

    After choosing a direction (from begin to end or vice versa, assumed to be from begin to end in
    the following), the clusters are iterated over in this direction, comparing always two
    neighbouring clusters, for simplicity labelled with 1 and 2. In each comparison the following 
    steps are performed
    - If it turns out, through TestClusters(Cluster *,Cluster *,int),
      that a shuffled cluster configuration is to be chosen the particles between the clusters
      are interchanged in the following way
      \f[
      C^{(1)}(q_1,\,\bar q_2)+C^{(2)}(q_3,\,\bar q_4) \longrightarrow
      C^{(1)}(q_1,\,\bar q_4)+C^{(2)}(q^{(2)}_3,\,\bar q_2)\,.
      \f]
    - Independent of whether the interchange took place or not, the two clusters' kinematical
      configurations/weights are compared, i.e. 
      \f[
      {\cal W} = \left\{\begin{array}{ll}
                 \frac{w_{14}}{w_{14}+w_{23}} & \mbox{\rm if shuffled} \\[2mm]
		 \frac{w_{12}}{w_{12}+w_{34}} & \mbox{\rm else} 
		 \end{array}\right.\,.
      \f]
      According to this weight, either cluster 1 or cluster 2 is then compared with the next cluster 
      in the sequence.
  */
  /*!
    \fn void Colour_Reconnections::Two_Singlet_CR(Cluster_List * cl1,Cluster_List * cl2)
    \brief Models colour reconnections between two singlet lists.

    \todo Not implemented yet.
  */
  /*!
    \fn bool Colour_Reconnections::TestClusters(Cluster * cl1,Cluster * cl2,int gen=1)
    \brief Returns "true" if two clusters are to be reshuffled, or "false" in the other case.

    This decision is based on a combined weight, consisting of the 
    KinematicWeight(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,const ATOOLS::Vec4D &) and the 
    ColourWeight(int).
  */
  /*!
    \fn double Colour_Reconnections::KinematicWeight(const ATOOLS::Vec4D & mom1,
                                                     const ATOOLS::Vec4D & mom2,
			                             const ATOOLS::Vec4D & mom3,
						     const ATOOLS::Vec4D & mom4)
    \brief The kinematic weight to decide whether two clusters are to be shuffled. It is based on some
    notion of being "kinematically favoured", i.e. closer in phase space.

    Labelling the four momenta by 1 to 4, where, by default, the pairs {1, 2} and {3, 4} belong to the two
    respctive clusters, the kinematical weight is given by the individual weights as
    \f[
    {\cal W}_k = \displaystyle{\frac{\frac{t_0}{t_0+(w_{14}+w_{23})^2}}
                                    {\frac{t_0}{t_0+(w_{12}+w_{34})^2}+\frac{t_0}{t_0+(w_{14}+w_{23})^2}}}
    \f] 
    If m_kinmode equals 1, the individual cluster weights are given by
    \f[
    w_{ij} = \sqrt{(p_i+p_j)^2}\,.
    \f]

    \todo Implement more cases, especially based on relative transverse momentum.
  */
  /*!
    \fn double Colour_Reconnections::ColourWeight(int gen=1)
    \brief The colour weight to decide whether two clusters are to be shuffled. So far, there is a trivial choice
    only.

    If m_colmode or the argument equal zero, the colour weight is zero, hence there is no reshuffling of clusters.
    If m_colmode equals 1, the colour weight is given by
    \f[
    {\cal W}_c = \frac19\,.
    \f]
  */

}
#endif
