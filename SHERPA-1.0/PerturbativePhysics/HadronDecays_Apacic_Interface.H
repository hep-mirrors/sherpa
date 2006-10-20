#ifndef HadronDecays_Apacic_Interface_H
#define HadronDecays_Apacic_Interface_H

#include "Perturbative_Interface.H"
#include "Tree.H"
#include "Blob.H"
#include "Poincare.H"
#include <list>

namespace SHERPA {

  class Interface_Tools;

  class HadronDecays_Apacic_Interface: public Perturbative_Interface {
  private:
    ATOOLS::Blob * p_blob;

    std::list<std::list<ATOOLS::Particle * > * > m_singlets;
    std::list<APACIC::Tree * >                   m_trees;

    ATOOLS::Poincare m_boost;               

    bool FillTree(std::list<ATOOLS::Particle *> *);
  public :

    // constructor
    HadronDecays_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);

    // destructor
    ~HadronDecays_Apacic_Interface();

    // member functions
    ATOOLS::Return_Value::code DefineInitialConditions(ATOOLS::Blob *blob);
    bool FillBlobs(ATOOLS::Blob_List *blobs); 

    int  PerformShowers();
    int  PerformDecayShowers();
  };
}
#endif
