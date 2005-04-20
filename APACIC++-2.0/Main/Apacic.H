#ifndef Apacic_H
#define Apacic_H

#include "Model_Base.H"
#include "Particle_List.H"
#include "Blob_List.H"
#include "Data_Read.H"

namespace ATOOLS { class Blob_Data_Base; }
namespace PDF    { class ISR_Handler;    }

namespace APACIC {

  class Initial_State_Shower;
  class Final_State_Shower;
  class Tree;

  //! The Interface between Matrix element and Parton Shower
  class Apacic {
    bool                   m_isron,m_fsron,m_showers;
    Initial_State_Shower * p_inishower;
    Final_State_Shower   * p_finshower;
    Tree                ** p_initrees;
    Tree                 * p_fintree;
    std::vector<double>    m_info_cms, m_info_lab;
  public :
    Apacic(PDF::ISR_Handler * isr,MODEL::Model_Base * model,int maxjetnumber,
	   bool isron,bool fsron,ATOOLS::Data_Read * dataread);
    ~Apacic(); 
    void PrepareTrees();
    void SetJetvetoPt2(const double q2i, const double q2f);
    void SetFactorisationScale(const double scale);
    int  PerformShowers(bool,bool,int,int,double,double,double ycut);
    bool ExtractPartons(bool,bool,ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);

    Tree  * FinTree();
    Tree ** IniTrees();

    Final_State_Shower   * FinShower();
    Initial_State_Shower * IniShower();

    ATOOLS::Blob_Data_Base *const Info(const int frame) const;


  };

  inline Tree  *Apacic::FinTree()  { return p_fintree; }
  inline Tree **Apacic::IniTrees() { return p_initrees; }

  inline Final_State_Shower   *Apacic::FinShower() { return p_finshower; }
  inline Initial_State_Shower *Apacic::IniShower() { return p_inishower; }

  /*!
    \file
    \brief contains the interface class APACIC::Apacic
  */

  // example see Run_Parameter.H
  /*!
    \mainpage
  */
}
#endif







