#include "CXXFLAGS.H"

#ifndef Ahadic_H
#define Ahadic_H


#include "Hadronisation_Parameters.H"
#include "Cluster_Formation_Handler.H"
#include "Cluster_Decay_Handler.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Particle_List.H"
#include "Return_Value.H"


namespace AHADIC {
  class Ahadic {
  private:
    bool                        m_fullinfo;
    int                         m_maxtrials;
    Cluster_Formation_Handler * p_cformhandler;
    Cluster_Decay_Handler     * p_cdechandler;
  public:
    Ahadic(std::string path=std::string("./"),
	   std::string file=std::string("Cluster.dat"),bool=false);
    ~Ahadic();

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob_List *);
  };
  /*!
    \file Ahadic.H
    \brief Contains the class AHADIC::Ahadic.
  */

  /*! 
    \class Ahadic
    \brief The main steering class of cluster fragmentation in SHERPA.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Ahadic
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Cluster_Formation_Handler * Ahadic::p_cformhandler
    \brief A pointer to the Cluster_Formation_Handler.
  */
  /*!
    \var Cluster_Decay_Handler * Ahadic::p_cdechandler
    \brief A pointer to the Cluster_Decay_Handler.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Ahadic
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Ahadic::Ahadic(std::string path=std::string("./"),
	               std::string file=std::string("Cluster.dat"),bool=false)
    \brief Initialises the main objects responsible for the fragmentation.

    In particular, these are the
    - Hadronisation_Parameters, <br>
      which contain also the Constituents, the All_Hadron_Multiplets,
      the All_Single_Transitions, and the All_Double_Transitions objects as well as the service class
      Momenta_Stretcher. For each of them, only one instance will be initialised in AHADIC.
    - Cluster_Formation_Handler, <br>
      which contains also the Gluon_Decayer, the Cluster_Former, the Colour_Reconnections and the
      Cluster_Transformer object. The pointer to the latter one will also be put as argument into the 
      Cluster_Decay_Handler.
    - Cluster_Decay_Handler, <br>
      which contains a Cluster_Decayer_Base, and through it the strategies for cluster decays as implemented in
      Cluster_Part and Hadron_Part.
  */
  /*!
    \fn Ahadic::~Ahadic()
    \brief Deletes the Cluster_Formation_Handler and the Cluster_Decay_Handler.
  */
  /*!
    \fn void Ahadic::Hadronize(ATOOLS::Blob_List *)
    \brief Forms and decays the clusters, filling in a fragmentation blob. 

    The basic methods for this are
    - Cluster_Formation_Handler::FormClusters(ATOOLS::Blob_List *) and
    - Cluster_Decay_Handler::DecayClusters(Cluster_List *,ATOOLS::Blob *).
  */
}
#endif
