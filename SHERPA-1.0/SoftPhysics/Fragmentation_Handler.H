#ifndef Fragmentation_Handler_h
#define Fragmentation_Handler_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Lund_Fortran_Interface.H"
#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Fragmentation_Handler {
    std::string              m_dir,m_file;
    std::string              m_fragmentationmodel;
    double                   m_lund_a,m_lund_b,m_lund_sigma;
    int                      m_mode;
    APHYTOOLS::Parton_List * p_mypartons;
    APHYTOOLS::Blob_List   * p_myblobs;
    APHYTOOLS::Blob        * p_myblob;
    Lund_Fortran_Interface * p_lund; 
    bool   ExtractSinglets(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    bool   FindConnected(APHYTOOLS::Blob_List *,APHYTOOLS::Parton *,
			 APHYTOOLS::Blob *);
  public :
    Fragmentation_Handler(std::string,std::string);
    ~Fragmentation_Handler();
    bool PerformFragmentation(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List * = NULL);
    void EmptyMyLists();

    Lund_Fortran_Interface * GetLundFortranInterface();
    std::string FragmentationModel() { return m_fragmentationmodel; }
  };
  /*!
    Interface to the Lund String. More fine-tuning needed.
  */
  /*!
    Used for recursive construction of singlet chains.
    Method compare colour component no 1 of parton with colour component no 2
    of list to fill singlets.
  */
}
#endif

