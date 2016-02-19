#ifndef AHADIC_Main_Ahadic_H
#define AHADIC_Main_Ahadic_H
#include "ATOOLS/Org/CXXFLAGS.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace AHADIC {


  class Cluster_Formation_Handler;
  class Cluster_Decay_Handler;

  class Ahadic {
  private:
    bool                        m_fullinfo;
    int                         m_maxtrials;
    std::string                 m_anadir;
    Cluster_Formation_Handler * p_cformhandler;
    Cluster_Decay_Handler     * p_cdechandler;
 
    Cluster_List m_clulist;
    Cluster_List m_prilist;

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob *,int retry=0);
    void Reset();
    bool SanityCheck(ATOOLS::Blob * blob,double norm=-1.);
    void CleanUp(ATOOLS::Blob * blob=NULL);
  public:
    Ahadic(std::string path=std::string("./"),
	   std::string file=std::string("Cluster.dat"));
    ~Ahadic();

    void writeHistos();
    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob_List *);

  };
  /*!
    \namespace AHADIC
    \brief SHERPA's version of the cluster fragmentation.
  */
  
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
	               std::string file=std::string("Cluster.dat"),bool ana=false)
    \brief Initialises the main objects responsible for the fragmentation.

    In particular, these are the
    - Hadronisation_Parameters, <br>
      which contain also the Constituents, the All_Hadron_Multiplets,
      the All_Single_Transitions, and the All_Double_Transitions objects as well as the service class
      Momenta_Stretcher. For each of them, only one instance will be initialised in AHADIC.
    - Cluster_Formation_Handler, <br>
      which contains also the Gluon_Decayer, the Cluster_Former, the and Colour_Reconnections.
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
