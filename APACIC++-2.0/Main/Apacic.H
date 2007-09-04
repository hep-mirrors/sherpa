#ifndef Apacic_H
#define Apacic_H

#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Tree.H"

#include "ISR_Handler.H"

#include "Model_Base.H"
#include "Particle_List.H"
#include "Blob_List.H"
#include "Data_Read.H"
#include "Jet_Veto.H"

namespace APACIC {

  class Apacic {
  private:

    Initial_State_Shower *p_inishower;
    Final_State_Shower   *p_finshower;

    Jet_Veto *p_jetveto;

    Tree **p_initrees;
    Tree  *p_fintree;

    ATOOLS::Vec4D m_cms;

    bool m_isron, m_fsron, m_showers, m_last_ljv;

    jv::mode m_jvm;

  public :
    
    // constructor
    Apacic(PDF::ISR_Handler *const isr,MODEL::Model_Base *const model,
	   ATOOLS::Jet_Finder *const jf,ATOOLS::Data_Read *const dataread);

    // destructor
    ~Apacic(); 

    // member functions
    int  PerformShowers(const double &x1=1.0,const double &x2=1.0);
    bool ExtractPartons(const bool ini,const bool fin,
			ATOOLS::Blob_List *const blist,
			ATOOLS::Particle_List *const plist=NULL);

    void PrepareTrees();
    void OutputTrees();

    // inline functions
    inline void SetFactorisationScale(const double &scale)
    { if (m_isron) p_inishower->SetFactorisationScale(scale); }

    inline Tree  *FinTree() const  { return p_fintree;  }
    inline Tree **IniTrees() const { return p_initrees; }

    inline Final_State_Shower   *FinShower() const { return p_finshower; }
    inline Initial_State_Shower *IniShower() const { return p_inishower; }

  };// end of class Apacic

}// end of namespace APACIC

#endif







