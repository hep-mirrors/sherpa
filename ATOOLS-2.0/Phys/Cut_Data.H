#ifndef Cut_Data_h
#define Cut_Data_h

#include "Flavour.H"

namespace ATOOLS {

  class Cut_Data {
    Flavour * fl;
  public:
    double ** cosmin;
    double ** cosmax;
    double ** scut;
    double ** scut_save;
    double *  energymin;
    double *  energymax;
    int       ncut;

    Cut_Data();
    ~Cut_Data();
    void Init(int,Flavour *);
    void Init(Cut_Data *,Flavour *);
    void Update(double,double);
  };    
}

#endif







