#ifndef Amegic_Apacic_Interface_H
#define Amegic_Apacic_Interface_H

#include "Perturbative_Interface.H"
#include "Cluster_Partons_Base.H"
#include "Tree_Filler.H"
#include "XS_Group.H" 
#include "XS_Base.H" 
#include "Blob.H"


namespace SHERPA {
  class Amegic_Apacic_Interface : public Perturbative_Interface {
  private :
    bool   m_isdecay, m_ckkwon;
    int    m_type, m_lastshowerveto;
    double m_ycut,m_jetscale;
    double m_qmin_i, m_qmin_f, m_scale;

    int    m_nin, m_nout;

    Cluster_Partons_Base  * p_cluster;
    Tree_Filler           * p_filler;
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
    ATOOLS::Return_Value::code DefineInitialConditions(ATOOLS::Blob *);
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

    Starting with a \f$2 \to n \f$ process a \f$2 \to 2 \f$ process is determined
    employing Cluster_Partons::ClusterConfiguration().

    In addition "ME PS Interface" blobs are prepared, to befilled into the event 
    record (if the event was excepted). 
  */

  /*!
    \fn bool Amegic_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob);
    \brief Calculates Sudakov weights and -- if accepted -- initializes
           the partonshower and returns 1, or -- if rejected -- returns 0, or 3 .. to be documented!!!

    The steps are the following
     -# cluster the matrix element configuration to a \f$2 \to 2 \f$ process
        using ClusterConfiguration()
     -# determine the starting scale and colour connections of this core process, 
        with the help of  EXTRAXS library, or via Cluster_Partons::SetColours()
     -# evaluate the NLL Sudakov weight used for reweighting the ME kinematics
     -# if necessary perform a rejection according to the weight
     -# if accepted setup the shower evolution
    .
  */
  /*! 
    \fn     bool Amegic_Apacic_Interface::FillBlobs(ATOOLS::Blob_List * blobs)
    \brief  inserts the neccessary ME PS merging blobs into the event record
  */
  /*! 
    \fn     bool Amegic_Apacic_Interface::CleanBlobs(ATOOLS::Blob_List * blobs)
    \brief  deletes all ME PS merging blobs in the event record
  */
  /*!
    \fn int Amegic_Apacic_Interface::PerformShowers()
    \brief  calls the appropriate shower routine using the Shower_Hander
  */

}
#endif







