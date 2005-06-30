#ifndef Amegic_Apacic_Interface_H
#define Amegic_Apacic_Interface_H

#include "Perturbative_Interface.H"

#include "Cluster_Partons.H"

#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Blob.H"

#include "Jet_Finder.H"


namespace SHERPA {
  class Amegic_Apacic_Interface : public Perturbative_Interface {
  private :
    bool   m_isdecay;
    int    m_type, m_lastshowerveto;
    double m_ycut,m_jetscale;
    double m_qmin_i, m_qmin_f, m_scale;

    int    m_nin, m_nout;

    ATOOLS::Jet_Finder * p_jf;
    Cluster_Partons    * p_cluster;
    /*! 
      Library of 2 -> 2 processes from Extra distribution (to extract the
      initial colour structure of the hardest 2->2 process within the hard blob).
      This library extends itself.
    */
    EXTRAXS::XS_Base      * p_xs; 
    EXTRAXS::XS_Group     * p_two2two, * p_one2N;

    ATOOLS::Blob * p_blob_psme_IS,* p_blob_psme_FS;

    bool ClusterConfiguration(ATOOLS::Blob * blob);
  public :
    Amegic_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~Amegic_Apacic_Interface(); 
    int  DefineInitialConditions(ATOOLS::Blob *);
    bool FillBlobs(ATOOLS::Blob_List *); 

    int  PerformShowers();
    int  PerformDecayShowers();
  };

  /*!
    \file 
    \brief Contains the class SHERPA::Amegic_Apacic_Interface
  */

  /*!
    \class Amegic_Apacic_Interface
    \brief Interface of ME and PS, where MEs are supplied by AMEGIC++ and
           PS is perform via APACIC++

    This is the central interface class, steering the various steps of the
    merging procedure. It is derived from the abstract class 
    Perturbative_Interface.
    The Amegic_Apacic_Interface is used through consecutive calls
    to the following methods:
    -# DefineInitialConditions(),
    -# ClusterConfiguration(),
    -# PerformShowers(),
    -# FillBlobs().

    Apart from ClusterConfiguration(), which is obsolete for
    instance for $2\to 2$ processes, these general methods have to be
    provided by any realisation of a Perturbative_-Interface.
  */
  /*!
    \var int m_type
    \brief Type of the process considered, important for the Jetfinder p_jf.

    -# m_type = 1: Lepton-Lepton collisions
    -# m_type = 4: Hadron-Hadron collisions (so far: anything else!)
  */

  /*!
    \fn  Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,Shower_Handler * shower);
    \brief Constructor, initializing ATOOLS::Jet_Finder, Cluster_Partons, EXTRAXS::XS_Group
  */

  /*!
    \fn Amegic_Apacic_Interface::~Amegic_Apacic_Interface(); 
    \brief Destructor
  */

  /*!
    \fn bool Amegic_Apacic_Interface::ClusterConfiguration(ATOOLS::Blob * blob);
    \brief Wrapper around the class Cluster_Partons, using a jet algoritm to 
           obtain a shower history.
	   
    This method is used to obtain a pseudo parton shower history, i.e. starting with a 
    \f$2 \to n \f$ process a \f$2 \to 2 \f$ process is determined. The clustering
    actually is achieved in the helper class, through the method 
    Cluster_Partons::ClusterConfiguration(). This method
    merely forms an intelligent wrapper around it, and it prepares
    merging Blobs of the type "ME PS Interface". These Blobs
    are used to translate the on-shell partons from the matrix element
    into the off-shell partons experiencing the parton shower. Therefore
    they are filled after the parton shower evolution.
  */

  /*!
    \fn bool Amegic_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob);
    \brief Calculates Sudakov weights and -- if accepted -- initializes
           the partonshower.
	   
    This method performs all steps necessary in order to construct a
    pseudo parton shower history and its corresponding weight, to 
    accept or reject this configuration, and, eventually, to initialise
    the parton shower. In particular for the first task, it heavily
    relies on the helper class Cluster_Partons. Depending on the success of the 
    procedure, the integer return value of this method is ``0'', ``1'' or ``3'' 
    indicating a rejected event, an accepted event, or a rejected event after the 
    lose-jet-veto, respectively. 

    The steps are the following:
     -# Cluster the matrix element configuration to a \f$2 \to 2 \f$ process
        using ClusterConfiguration().
     -# Determine the starting scale and colour connections of this core process, 
        with the help of EXTRAXS library, or via Cluster_Partons::SetColours().
     -# evaluate the NLL Sudakov weight through Cluster_Partons::CalculateWeight()}.
        It is used for reweighting the ME kinematics.
     -# If accepted initialise the parton shower evolution by employing 
        Cluster_Partons::FillTrees().
  */
  /*! 
    \fn     bool Amegic_Apacic_Interface::FillBlobs(ATOOLS::Blob_List * blobs)
    \brief  inserts the neccessary ME PS merging blobs into the event record
  */
  /*! 
    \fn     bool Amegic_Apacic_Interface::CleanBlobs(ATOOLS::Blob_List * blobs)
    \brief  deletes all ME PS merging blobs in the event record -- obsolete
  */
  /*!
    \fn int Amegic_Apacic_Interface::PerformShowers()
    \brief  calls the appropriate shower routine using the Shower_Hander

    This method calls the appropriate routine in the Shower_Hander. 
    The jet-veto scale and the renormalisation scale, which have
    been determined in the merging procedure before, are
    handed over to the parton shower, also through the 
    Shower_Hander (Remember, they may change because of, e.g., the highest
    multiplicity treatment described elsewhere.). If the shower evolution was successful,
    FillBlobs() inserts the prepared and filled ``ME PS Interface'' and the shower
    Blobs into the event record.
   */
}
#endif







