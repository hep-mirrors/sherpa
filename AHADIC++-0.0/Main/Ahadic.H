#ifndef Ahadic_H
#define Ahadic_H

#include "Cluster.H"
#include "Clusters_2_Hadrons.H"
#include "Cluster_Decay_Handler.H"
#include "Cluster_Formation_Handler.H"

#include "Particle_List.H"
#include "Blob_List.H"


namespace AHADIC {
  class Ahadic {
  private:
    std::string                 m_path, m_file;

    Cluster_List              * p_clusters;
    ConstituentMasses         * p_constituents;
    Hadron_List               * p_hadrons;
    Clusters_2_Hadrons        * p_c2hadrons;
    Cluster_Formation_Handler * p_cformer;
    Cluster_Decay_Handler     * p_cdecayer;
  public:
    Ahadic(std::string _m_path,std::string _m_file);
    bool Hadronize(ATOOLS::Blob_List *bloblist,
		   ATOOLS::Particle_List *pl=NULL);
  };
}
#endif
