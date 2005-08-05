#ifndef SimpleXS_CSS_Interface_h
#define SimpleXS_CSS_Interface_h

#include "Perturbative_Interface.H"

namespace SHERPA {
  class SimpleXS_CSS_Interface : public Perturbative_Interface {
  private:
    ATOOLS::Blob * p_blob;
    bool InitColours(ATOOLS::Blob *);
    bool Valid(ATOOLS::Blob *);
  public :
    // constructor
    SimpleXS_CSS_Interface(Matrix_Element_Handler *,Shower_Handler *);
    // destructor
    ~SimpleXS_CSS_Interface();

    // member functions
    int  DefineInitialConditions(ATOOLS::Blob *blob);
    bool FillBlobs(ATOOLS::Blob_List *blobs); 

    int  PerformShowers();
    int  PerformDecayShowers() { return 1; }
  };
}
#endif
