#ifndef SimpleXS_CSS_Interface_h
#define SimpleXS_CSS_Interface_h

#include "Perturbative_Interface.H"
#include "Poincare.H"



namespace SHERPA {
  class SimpleXS_CSS_Interface : public Perturbative_Interface {
  protected :
    int m_hadron;
    int m_lepton;
  private:
    ATOOLS::Vec4D   *p_momenta;
    ATOOLS::Flavour *p_flavours;
    
    ATOOLS::Blob * p_blob, * p_hard;
    ATOOLS::Poincare m_boost;
    ATOOLS::Poincare m_zrot;
  
    double m_x, m_scale;
    
    bool InitColours(ATOOLS::Blob *);
    bool InitBreitFrame(ATOOLS::Blob *);
    void BoostInBreitFrame(ATOOLS::Vec4D&);
    bool Valid(ATOOLS::Blob *);
    ATOOLS::Blob* ConstructSinglets(ATOOLS::Blob *,const double);
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
