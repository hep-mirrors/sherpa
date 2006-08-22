#ifndef Cluster_Transformer_H
#define Cluster_Transformer_H

#include "Transitions.H"
#include "Cluster.H"
#include "Blob.H"
#include "Return_Value.H"

namespace AHADIC {
  struct ctrans {
    enum code {
      nothing        = 0,
      forceddecay    = 1,
      photonemission = 10,
      pi0emission    = 11
    };
  };

  class Cluster_Transformer {
  private:
    ctrans::code             m_mode;
    double                   m_offset;
    double                   m_gammaenergy;
    All_Single_Transitions * p_stransitions;
    All_Double_Transitions * p_dtransitions;

    void DecayCluster(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,ATOOLS::Blob *);
  public:
    Cluster_Transformer();
    ~Cluster_Transformer();
    ATOOLS::Return_Value::code TreatSingleCluster(Cluster *,ATOOLS::Blob *);
    ATOOLS::Return_Value::code TreatClusterList(Cluster_List *,ATOOLS::Blob *);
  };

  /*!
    \file 
    \brief Contains the class AHADIC::Cluster_Transformer and the corresponding struct AHADIC::ctrans 
  */

  /*!
    \struct ctrans
    \brief Short description of cluster transformation type, used to steer the operational mode
    of AHADIC::Cluster_Transformer. Available labels so far are: 
    nothing (0), forceddecay (1), pi0emission (2).
  */
  /*! 
    \class Cluster_Transformer
    \brief This class deals with single light clusters.

    There are different modes for their decays, available/foreseen so far are:
    - nothing: Nothing happens to those clusters which must be transformed into a hadron.
      This mode is implemented for testing purposes only.
    - forceddecay: In this mode, a AHADIC::Double_Transition_Map will be used for isotropic decays
      including the pseudoscalar octet. This still waits for implementation.
    - pi0emission: Under the emission of a \f$ \pi_0\f$ the clusters transform into a hadron 
      with the same flavour quantum numbers.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Transformer
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var ctrans::code Cluster_Transformer::m_mode
    \brief The mode which determines how single light clusters are dealt with. 
  */
  /*!
    \var double Cluster_Transformer::m_offset
    \brief The minimal offset between a cluster and the heaviest hadron with the same flavour 
    composition before the cluster must be trasnformed into hadrons.
  */
  /*!
    \var All_Single_Transitions * Cluster_Transformer::p_transitions
    \brief A pointer to the AHADIC::All_Single_Transitions class. With tihs class checks are performed whether
    a cluster must decay into hadrons and a suitable hadron is found.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Transformer
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Transformer::Cluster_Transformer()
    \brief Default constructor. By default, so far, the transformation mode is set to pi0emission,
    the offset parameter and the AHADIC::All_Single_Transitions -pointer are taken from the 
    AHADIC::Hadronisation_Parameters.
  */
  /*!
    \fn Cluster_Transformer::~Cluster_Transformer()
    \brief The default destructor, so far an empty envelope.
  */
  /*!
    \fn bool Cluster_Transformer::TreatSingleCluster(Cluster_List *,ATOOLS::Blob *)
    \brief Deals with the first cluster of a list; called from the AHADIC::Cluster_Formation_Handler. By construction
    so far this method is applied on lists with one member only. Returns true only if a decay of the cluster
    into hadrons must be performed, otherwise it returns false. 

    If the first cluster of a list must and can be decayed into hadrons, as checked by the method
    AHADIC::All_Single_Transitions::MustTransit(Cluster *,ATOOLS::Flavour &,const double=0.), the hadrons are added as 
    outgoing particles to the Blob, the cluster is erased from the list and true is returned. True is also 
    returned, if the cluster must decay into hadrons but the operation mode is not to perform any decay (m_mode=no).
    In such a case, the cluster in question will be added to the shortest cluster singlet list in the
    AHADIC::Cluster_Formation_Handler, to allow for suitable four-momentum transfers. If the cluster does not need
    to be decayed immediately into hadrons, false is returned.
  */
  /*!
    \fn bool Cluster_Transformer::TreatSingleCluster(Cluster *,Part_List *)
    \brief Deals with single clusters; called from the AHADIC::Cluster_Decay_Handler. Returns true only if a decay 
    of the cluster into hadrons must be performed, otherwise it returns false. 
    
    The method is widely similar to Cluster_Transformer::TreatSingleCluster(Cluster_List *,ATOOLS::Blob *),
    the only difference being that emerging particles are added to a Part_List rather than to a Blob.
  */
  /*!
    \fn void Cluster_Transformer::TreatClusterList(Cluster_List *,ATOOLS::Blob *)
    \brief Deals with full cluster lists and checks whether any of the clusters must immediately transform into
    hadrons. If so, the transformations are performed, resulting usually in a shift of mass and momenta. In such
    a case the momenta of all remaining objects - clusters or hadrons - are shifted accordingly.

    This method is called from the Cluster_Formation_Handler, passing all primordial cluster lists through
    it. For each cluster it is checked, whether its mass is smaller than the mass of the heaviest hadron with
    same flavour composition plus the offset, Cluster_Transformer::m_offset. This check is performed through the method
    AHADIC::All_Single_Transitions::MustTransit(Cluster *,ATOOLS::Flavour &,const double=0.). If any transformations 
    occur, the method AHADIC::Hadronisation_Parameters::AdjustMomenta(const int,ATOOLS::Vec4D *,const double *) is 
    employed to shift the momenta of all remaining objects on the right mass-shell. Finally, all 
    emerging primoridal hadrons are added as outgoing particles to the blob.
  */
  /*!
    \fn void Cluster_Transformer::DecayCluster(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,ATOOLS::Blob *)
    \brief Performs an isotropic decay of the cluster into the two (hadron-) flavours and adds the new particles
    to the blob.
  */
  /*!
    \fn void Cluster_Transformer::DecayCluster(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,Part_List *)
    \brief Performs an isotropic decay of the cluster into the two (hadron-) flavours and adds the new particles
    to the particle list.
  */
}
#endif
