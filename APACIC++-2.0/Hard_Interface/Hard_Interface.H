#ifndef HardInterface_h
#define HardInterface_h

#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Interface_Tools.H"
#include "Tree.H"

#include "XS_Selector.H" 
#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Parton_List.H"
#include "Blob_List.H"

#include "ISR_Handler.H"

namespace APACIC {
  //! The Interface between Matrix element and Parton Shower
  class Hard_Interface {
    int nin,nout,maxjetnumber;

    ISR::ISR_Handler     * isr;
    /*! 
      Library of 2 -> 2 processes from Extra distribution (to extract the
      initial colour structure of the hardest 2->2 process within the hard blob).
      This library extends itself.
    */
    EXTRAXS::XS_Group    * two2two;
    EXTRAXS::XS_Selector * xsselector;

    /*!
      Pointers to the parton showers and the trees. 
      Note that for incoming parton showers you have two trees.
    */
    Initial_State_Shower * ini_shower;
    Final_State_Shower   * fin_shower;
    Tree                 * fin_tree;
    Tree                ** ini_trees;

    Interface_Tools      * tools;

    bool InitColours(APHYTOOLS::Blob *);
    bool InitialConditions(APHYTOOLS::Blob *,EXTRAXS::XS_Base *);
  public :
    Hard_Interface(ISR::ISR_Handler *&,int,bool); 
    ~Hard_Interface(); 
    void PrepareTrees();
    bool ReduceBlob(APHYTOOLS::Blob *);
    int  PerformShowers(bool,bool);
    bool ExtractPartons(bool,bool,APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    void FinalStats();

    Tree             * FinTree()  { return fin_tree; }
    Tree            ** IniTrees() { return ini_trees; }
    Interface_Tools  * Tools()    { return tools; } 
  };
}
#endif







