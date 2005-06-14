#ifndef Ahadic_H
#define Ahadic_H


#include "Hadronisation_Parameters.H"
#include "Cluster_Formation_Handler.H"
#include "Cluster_Decay_Handler.H"
#include "Cluster.H"
#include "Blob_List.H"


namespace AHADIC {
  class Ahadic {
  private:
    Cluster_Formation_Handler * p_cformhandler;
    Cluster_Decay_Handler     * p_cdechandler;
  public:
    Ahadic(std::string path=std::string("./"),
	   std::string file=std::string("Cluster.dat"),bool=false);
    ~Ahadic();

    void Hadronize(ATOOLS::Blob_List *);
  };
}
#endif
