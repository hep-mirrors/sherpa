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
    int            m_N_of_singlets;
    ATOOLS::Blob * p_blob;

    std::vector<ATOOLS::Particle * > m_particles;

    ATOOLS::Poincare m_boost;               

    bool FillTree(APACIC::Tree *);
    bool FillBinaryDecayTree(APACIC::Tree *);
    bool FillTertiaryDecayTree(APACIC::Tree *);
    bool FillSpectatorDecayTree(APACIC::Tree *);
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
