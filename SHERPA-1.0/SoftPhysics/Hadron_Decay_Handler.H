#ifndef Hadron_Decay_Handler_h
#define Hadron_Decay_Handler_h


#include <set> 
#include <iterator> 
#include <iostream> 

#include "Lund_Fortran_Interface.H"
#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Hadron_Decay_Handler {
  protected :
    std::string              m_dir,m_file;
    Lund_Fortran_Interface * p_lund; 
    bool ReconstructLundHadronDecays();
  public :
    Hadron_Decay_Handler(string,string);
    Hadron_Decay_Handler(string,string,Lund_Fortran_Interface *);
    ~Hadron_Decay_Handler();
    bool FillHadronDecayBlobs(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
  };
}
#endif

