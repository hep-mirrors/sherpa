#ifndef ME_PSInterface_h
#define ME_PSInterface_h

#include "Hard_Interface.H"

#include "Cluster_Partons.H"

#include "XS_Selector.H" 
#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Parton_List.H"
#include "Blob_List.H"
#include "Blob.H"

#include "Jet_Finder.H"
#include "ISR_Handler.H"

namespace SHERPA {
  //! The Interface between Matrix element and Parton Shower
  class ME_PS_Interface {
    bool   ini,fin;
    bool   NLLweight;
    double ycut,jetscale;
    int    nin,nout,maxjetnumber;

    APHYTOOLS::Flavour     * fl;
    AMATOOLS::Vec4D        * p;
    APHYTOOLS::Jet_Finder  * jf;
    ISR::ISR_Handler       * isr;
    Cluster_Partons        * cluster;
    /*! 
      Library of 2 -> 2 processes from Extra distribution (to extract the
      initial colour structure of the hardest 2->2 process within the hard blob).
      This library extends itself.
    */
    EXTRAXS::XS_Base       * xs; 
    EXTRAXS::XS_Group      * two2two;
    EXTRAXS::XS_Selector   * xsselector;


    APACIC::Hard_Interface * ps;
    APHYTOOLS::Blob        * blob;
  public :
    ME_PS_Interface(ISR::ISR_Handler *,int); 
    ~ME_PS_Interface(); 
    bool Treat(AMEGIC::Process_Base *,APHYTOOLS::Blob *,int);
    int  PerformShower(AMEGIC::Process_Base *,int);
    bool ExtractPartons(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    bool ExtractFinalPartons(APHYTOOLS::Parton_List *);
  };
}
#endif







