#ifndef Cluster_Former_H
#define Cluster_Former_H

#include "Cluster.H"

namespace AHADIC {
  class Cluster_Former {
  private:
  public:
    Cluster_Former();
    ~Cluster_Former();
    void ConstructClusters(SP(Proto_Particle_List),Cluster_List *);
  };
  /*!
    \file Cluster_Former.H
    \brief Contains the class AHADIC::Cluster_Former. 
  */

  /*! 
    \class Cluster_Former
    \brief Constructs a cluster list from a particle list, through the (so far) only relevant method
    void ConstructClusters(Part_List * plin,Cluster_List * clout).
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Former
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Former
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Former::Cluster_Former()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Cluster_Former::~Cluster_Former()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Cluster_Former::ConstructClusters(Part_List *,Cluster_List *)
    \brief Constructs a cluster list from a particle list.

    It iterates over the particles in the list and puts particles with colour configuration (a,0) into 
    clusters together with particles having colour configuration (0,a). In case the particle list is 
    empty or if there are unwanted colour configurations in the list, the program will abort.
  */

}
#endif
