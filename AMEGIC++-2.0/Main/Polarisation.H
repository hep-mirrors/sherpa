#ifndef Polarisation_h
#define Polarisation_h

#include "Basic_Sfuncs.H"
#include "Single_Amplitude.H"

namespace AMEGIC {
  //! Contains everything with massless and massive polarisations.
  class Polarisation {
    double Mass_Norm; 
    int nmass,npol;
    int** mass_pol;
    int no;
  public:
    Polarisation();
    ~Polarisation();
    double Spin_Average(int,APHYTOOLS::Flavour*);
    int    Massless_Vectors(int,APHYTOOLS::Flavour*);
    int    Massive_Vectors(int,APHYTOOLS::Flavour*);
    void   Add_Extern_Polarisations(Basic_Sfuncs*,APHYTOOLS::Flavour*,Helicity*);
    void   Attach(int,APHYTOOLS::Flavour*);
    void   Reset_Gauge_Vectors(int,AMATOOLS::vec4d*,AMATOOLS::vec4d);
    void   Set_Gauge_Vectors(int,AMATOOLS::vec4d*,AMATOOLS::vec4d);
    double Massless_Norm(int,APHYTOOLS::Flavour*,Basic_Sfuncs*);
    double Massive_Norm() { return Mass_Norm; }
    void   Replace_Numbers(int,APHYTOOLS::Flavour*,Single_Amplitude*);
  };
}
#endif




