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

namespace MOCAIC {
  class Soft_Interface {
    ATOOLS::Flavour ** constituents;
    int     *  n_const;

    /*!
      Intrinsic partonlist: After ExtractSinglets it will be colour-ordered.
    */
    ATOOLS::Parton_List * mypartons;
    ATOOLS::Blob_List   * myblobs;
    ATOOLS::Blob        * myblob;
    /*!
      Interface to the Lund String. More fine-tuning needed.
    */
    Lund_Fortran_Interface * lund; 
    /*!
      Used for recursive construction of singlet chains.
      Method compare colour component no 1 of parton with colour component no 2
      of list to fill singlets.
    */
    bool   ExtractSinglets(ATOOLS::Blob_List *,ATOOLS::Parton_List *);
    bool   FindConnected(ATOOLS::Parton_List *,ATOOLS::Parton *,
			 ATOOLS::Blob *);
    int    Constituents(ATOOLS::Flavour,ATOOLS::Flavour *); 
    bool   FillHadron(ATOOLS::Parton_List *,ATOOLS::Blob *,ATOOLS::Parton *);
    double GetX(ATOOLS::Flavour,ATOOLS::Flavour,double);
  public :
    Soft_Interface();
    ~Soft_Interface();
    /*!
      Generic methods to extract singlets from the parton/bloblist
    */
    bool PerformFragmentation(ATOOLS::Blob_List *,ATOOLS::Parton_List *);

    bool HadronsToPartons(ATOOLS::Blob_List *,ATOOLS::Parton_List *);
    // bool ExtractSinglets(ATOOLS::Blob_List *);
    void EmptyMyLists();
  };

  //----------------------------------------------------------------------
  //    inline functions
  //----------------------------------------------------------------------

}
#endif

