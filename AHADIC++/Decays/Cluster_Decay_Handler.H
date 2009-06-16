#ifndef AHADIC_Decays_Cluster_Decay_Handler_H
#define AHADIC_Decays_Cluster_Decay_Handler_H

#include "AHADIC++/Decays/Cluster_Part.H"
#include "AHADIC++/Tools/Dipole_Splitter.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Decays/Cluster_Decay_Analysis.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    Soft_Cluster_Handler   * p_softclusters;
    Cluster_Part           * p_clus;
    Cluster_List           * p_clulist;

    Cluster_Decay_Analysis * p_analysis;
 
    ATOOLS::Blob * DecayIt(Cluster *);
    ATOOLS::Blob * ClusterDecayBlob(Cluster *,Cluster_List *);
  public:
    Cluster_Decay_Handler(Cluster_List *,bool=false);
    ~Cluster_Decay_Handler();

    int DecayClusters(ATOOLS::Blob *);
  };
  /*!
    \file Cluster_Decay_Handler.H
    \brief Contains the class AHADIC::Cluster_Formeration_Handler.
  */
  /*!
    \class Cluster_Decay_Handler
    \brief The Cluster_Decay_Handler organises the decay of clusters at the end of
    the parton shower phase.  

    The central (and also the only public) method of this class is 
    Cluster_Decay_Handler::DecayClusters(Cluster_List *,ATOOLS::Blob_List *),
    which is being called, once the previous cluster formation from partons, 
    managed by the class Cluster_Formation_Handler, is finished.  
  */
  /*!
    \var Soft_Cluster_Handler * p_softclusters
    \brief A pointer to the common Soft_Cluster_Handler, used in both the formation
    and decay phase to decide on forced cluster to hadron transitions and decays
    and to perform them.
  */
  /*!
    \var Cluster_Part * p_clus
    \brief Pointer to the Cluster_Part, responsible for the definition of the actual
    cluster decay kinematics.
  */
  /*!
    \var Cluster_Decay_Analysis * p_analysis
    \brief Pointer to the Cluster_Decay_Analysis.  
  */
  /*!
    \fn ATOOLS::Blob * Cluster_Decay_Handler::DecayIt(Cluster *)
    \brief Managing for the decay of one single cluster.

    The following steps are taken:
    - a decay blob (status: blob_status::needs_hadrondecays) is initialised
      and the cluster is added as incoming particle;
    - the method Cluster_Part::TestDecay(Cluster*) decides upon a flavour to be 
      popped out of the vacuum and instantiates two new clusters with according
      flavour content.  They are then linked as left and right pointers to the 
      decaying one;
    - finally, the cluster decay blob (without the two outgoing clusters added)
      is returned.    
  */
  /*!
    \fn Cluster_Decay_Handler::Cluster_Decay_Handler(Soft_Cluster_Handler *,bool=false)
    \brief Constructor of the class.  Initialises a new Cluster_Part and, eventually,
    a Cluster_Decay_Analysis object. 

    The boolean steers, whether the Cluster_Decay_Analysis object is instantiated.
  */
  /*!
    \fn Cluster_Decay_Handler::~Cluster_Decay_Handler()
    \brief Deletes the Cluster_Part and the  Cluster_Decay_Analysis object.
  */
  /*!
    \fn int Cluster_Decay_Handler::DecayClusters(Cluster_List *,ATOOLS::Blob_List * = NULL)
    \brief Managing the decays of the clusters in a Cluster_List into secondary
    clusters and, ultimately, the primordial hadrons.

    The idea is the following: each decaying cluster is connected with a cluster list,
    filled by its offsprings.  These emerge from the call to the method 
    Cluster_Decay_Handler::DecayIt(Cluster *) and, at that point, are still clusters.
    Then the Soft_Cluster_Handler tests the two offspring clusters, whether they must
    decay or transit into hadrons.  This is done through the method
    Soft_Cluster_Handler::TreatClusterList(Cluster_List *,ATOOLS::Blob *).  Now the
    offspring list is iterated over.  Each of them is treated in the following way:
    - First, the ATOOLS::Particle related to it is added as outgoing particle to the decay 
      ATOOLS::Blob of the original cluster.  
    - If the offspring does not have any filled left and right pointers, it did not decay into 
      hadrons and hence it can safely be added to the original cluster list for further treatment.  
    - If, on the other hand, the offspring does have filled left and right pointers, it decayed 
      into hadrons or clusters.  Then the offsprings of the offspring are tested: If they are 
      not hadrons, they are added to the original cluster list for further treatment.  

    Finally, the original cluster is deleted and taken out of the original cluster list.  

    This procedure is repeated until the original cluster list is empty, i.e., until all
    clusters are decayed into hadrons.
  */
}
#endif
