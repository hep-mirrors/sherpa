#ifndef SimpleXS_Apacic_Interface_h
#define SimpleXS_Apacic_Interface_h

#include "Perturbative_Interface.H"
#include "Interface_Tools.H"

namespace SHERPA {
  class SimpleXS_Apacic_Interface : public Perturbative_Interface {
    Interface_Tools  * p_tools;

    bool InitColours(APHYTOOLS::Blob *);
  public :
    SimpleXS_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~SimpleXS_Apacic_Interface();
    bool   DefineInitialConditions(APHYTOOLS::Blob *);
  };
}
#endif







