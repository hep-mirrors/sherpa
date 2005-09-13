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
#include "Jet_Finder.H"

namespace APACIC {
  //! The Interface between Matrix element and Parton Shower
  class Apacic {
    bool                   m_isron,m_fsron,m_showers;
    Initial_State_Shower * p_inishower;
    Final_State_Shower   * p_finshower;
    Tree                ** p_initrees;
    Tree                 * p_fintree;
  public :
    Apacic(PDF::ISR_Handler * isr,MODEL::Model_Base * model,ATOOLS::Jet_Finder *jf,ATOOLS::Data_Read * dataread);
    ~Apacic(); 
    int  PerformShowers(int=0,int=0,double=1.,double=1.,double ycut=1.e-2);
    bool ExtractPartons(bool,bool,ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
    bool FinalFSJetCheck(const double);
    void BoostInCorrectSystem();

    inline void PrepareTrees();
    inline void SetJetvetoPt2(const double q2i, const double q2f);
    inline void SetFactorisationScale(const double scale);
    inline void OutputTrees();
    inline ATOOLS::Vec4D PrepareFSR();

    Tree                 * const FinTree()   const { return p_fintree;   }
    Tree                ** const IniTrees()  const { return p_initrees;  }
    Final_State_Shower   * const FinShower() const { return p_finshower; }
    Initial_State_Shower * const IniShower() const { return p_inishower; }

    inline void SetMaxJetNumber(const size_t maxjets) 
    { 
      p_finshower->SetMaxJetNumber(maxjets); 
    }
  
  };

  void Apacic::PrepareTrees() {
    if (m_fsron) { p_fintree->Reset(); } 
    if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
  }

  void Apacic::SetJetvetoPt2(const double q2i, const double q2f)
  { 
    if (m_fsron) p_finshower->SetJetvetoPt2(q2f); 
    if (m_isron) p_inishower->SetJetvetoPt2(q2i,q2i); 
  }

  void Apacic::SetFactorisationScale(const double scale)
  {
    if (m_isron) p_inishower->SetFactorisationScale(scale);
  }

  void Apacic::OutputTrees() {
    if (m_fsron) 
      p_finshower->OutputTree(p_fintree);
    if (m_isron) {
      p_inishower->OutputTree(p_initrees[0]);
      p_inishower->OutputTree(p_initrees[1]);
    }
  }

  ATOOLS::Vec4D Apacic::PrepareFSR() {
    ATOOLS::Vec4D cms(p_fintree->GetRoot()->part->Momentum());
    ATOOLS::Poincare cmsb(cms);
    p_fintree->BoRo(cmsb);
    return cms;
  }
}
#endif







