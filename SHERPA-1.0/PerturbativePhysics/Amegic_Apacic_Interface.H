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
    double m_type;
    double m_ycut,m_jetscale;

    APHYTOOLS::Jet_Finder * p_jf;
    Cluster_Partons       * p_cluster;
    /*! 
      Library of 2 -> 2 processes from Extra distribution (to extract the
      initial colour structure of the hardest 2->2 process within the hard blob).
      This library extends itself.
    */
    EXTRAXS::XS_Base      * p_xs; 
    EXTRAXS::XS_Group     * p_two2two;

    APHYTOOLS::Blob * p_blob_psme;

    bool ClusterConfiguration(APHYTOOLS::Blob * blob);
  public :
    Amegic_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~Amegic_Apacic_Interface(); 
    bool   DefineInitialConditions(APHYTOOLS::Blob *);
    bool   FillBlobs(APHYTOOLS::Blob_List *); 
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
    \fn  Amegic_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    \brief Constructor 
  */

  /*!
    \fn ~Amegic_Apacic_Interface(); 
    \brief Destructor
  */

  /*!
    \fn bool Amegic_Apacic_Interface::ClusterConfiguration(APHYTOOLS::Blob * blob);
    \brief Wrapper around the class Cluster_Partons, using a jet algoritm to 
           obtain a shower history.
  */

  /*!
    \fn bool Amegic_Apacic_Interface::DefineInitialConditions(APHYTOOLS::Blob *);
    \brief Calculates Sudakov weights and -- if accepted -- initializes
           the partonshower and returns 1, or -- if rejected -- returns 0.
  */

}
#endif







