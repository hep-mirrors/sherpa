#ifndef Jet_Evolution_h
#define Jet_Evolution_h

#include "Event_Phase_Handler.H"
#include "Matrix_Element_Handler.H"
#include "Shower_Handler.H"
#include "Perturbative_Interface.H"

namespace SHERPA {
  class Jet_Evolution : public Event_Phase_Handler {
  private :
    Matrix_Element_Handler * p_mehandler;
    Shower_Handler         * p_showerhandler;
    Perturbative_Interface * p_interface;
  public :
    Jet_Evolution(Matrix_Element_Handler *,Shower_Handler *);
    ~Jet_Evolution();
    bool Treat(APHYTOOLS::Blob_List *);
    void CleanUp();
  };
}
#endif

