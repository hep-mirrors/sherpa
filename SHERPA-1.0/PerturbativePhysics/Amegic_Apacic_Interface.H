#ifndef Amegic_Apacic_Interface_h
#define Amegic_Apacic_Interface_h

#include "Perturbative_Interface.H"

//#include "Hard_Interface.H"
#include "Cluster_Partons.H"

#include "XS_Selector.H" 
#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Parton_List.H"
#include "Blob_List.H"
#include "Blob.H"

#include "XS_Selector.H" 
#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Jet_Finder.H"
#include "Blob.H"


namespace SHERPA {
  class Amegic_Apacic_Interface : public Perturbative_Interface {
  private :
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
    EXTRAXS::XS_Selector  * p_xsselector;

    APHYTOOLS::Blob       * p_blob;

    bool ClusterConfiguration(APHYTOOLS::Blob * blob);
  public :
    Amegic_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~Amegic_Apacic_Interface(); 
    bool   DefineInitialConditions(APHYTOOLS::Blob *);
  };
}
#endif







