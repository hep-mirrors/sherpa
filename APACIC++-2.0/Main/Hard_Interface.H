#ifndef HardInterface_h
#define HardInterface_h

#include "Model_Base.H"

#include "Parton_List.H"
#include "Blob_List.H"

#include "Data_Read.H"


namespace PDF {
  class ISR_Handler;
}

namespace APACIC {
  class Initial_State_Shower;
  class Final_State_Shower;
  class Tree;

  //! The Interface between Matrix element and Parton Shower
  class Hard_Interface {
    bool                   m_isron,m_fsron,m_showers;
    Initial_State_Shower * p_inishower;
    Final_State_Shower   * p_finshower;
    Tree                 * p_fintree;
    Tree                ** p_initrees;
  public :
    Hard_Interface(PDF::ISR_Handler *,MODEL::Model_Base *,int,bool,bool,AORGTOOLS::Data_Read *); 
    ~Hard_Interface(); 
    void PrepareTrees();
    int  PerformShowers(bool,bool,bool);
    bool ExtractPartons(bool,bool,APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List * = NULL);

    Tree             * FinTree()  { return p_fintree; }
    Tree            ** IniTrees() { return p_initrees; }
  };
}
#endif







