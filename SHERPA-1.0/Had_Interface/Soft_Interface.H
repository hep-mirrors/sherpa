#ifndef SoftInterface_h
#define SoftInterface_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Lund_Fortran_Interface.H"
#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Soft_Interface {
    APHYTOOLS::Flavour ** constituents;
    int     *  n_const;

    /*!
      Intrinsic partonlist: After ExtractSinglets it will be colour-ordered.
    */
    APHYTOOLS::Parton_List * mypartons;
    APHYTOOLS::Blob_List   * myblobs;
    APHYTOOLS::Blob        * myblob;
    /*!
      Interface to the Lund String. More fine-tuning needed.
    */
    Lund_Fortran_Interface * lund; 
    /*!
      Used for recursive construction of singlet chains.
      Method compare colour component no 1 of parton with colour component no 2
      of list to fill singlets.
    */
    bool   ExtractSinglets(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    bool   FindConnected(APHYTOOLS::Parton_List *,APHYTOOLS::Parton *,
			 APHYTOOLS::Blob *);
    int    Constituents(APHYTOOLS::Flavour,APHYTOOLS::Flavour *); 
    bool   FillHadron(APHYTOOLS::Parton_List *,APHYTOOLS::Blob *,APHYTOOLS::Parton *);
    double GetX(APHYTOOLS::Flavour,APHYTOOLS::Flavour,double);
  public :
    Soft_Interface();
    ~Soft_Interface();
    /*!
      Generic methods to extract singlets from the parton/bloblist
    */
    bool PerformFragmentation(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);

    bool HadronsToPartons(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    // bool ExtractSinglets(APHYTOOLS::Blob_List *);
    void EmptyMyLists();
  };

  //----------------------------------------------------------------------
  //    inline functions
  //----------------------------------------------------------------------

}
#endif

