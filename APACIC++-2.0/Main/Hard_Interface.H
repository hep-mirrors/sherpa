#ifndef HardInterface_h
#define HardInterface_h

#include "Model_Base.H"

#include "Particle_List.H"
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
    Tree                ** p_initrees;
    Tree                 * p_fintree;
  public :
    Hard_Interface(PDF::ISR_Handler *,MODEL::Model_Base *,int,bool,bool,ATOOLS::Data_Read *); 
    ~Hard_Interface(); 
    void PrepareTrees();
    void SetJetvetoPt2(const double pt2);
    int  PerformShowers(bool,bool,int,double,double);
    bool ExtractPartons(bool,bool,ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);

    Tree  * FinTree();
    Tree ** IniTrees();

    Final_State_Shower   * FinShower();
    Initial_State_Shower * IniShower();

  };

  inline Tree  *Hard_Interface::FinTree()  { return p_fintree; }
  inline Tree **Hard_Interface::IniTrees() { return p_initrees; }

  inline Final_State_Shower   *Hard_Interface::FinShower() { return p_finshower; }
  inline Initial_State_Shower *Hard_Interface::IniShower() { return p_inishower; }

}
#endif







