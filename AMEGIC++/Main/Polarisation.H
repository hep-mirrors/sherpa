#ifndef AMEGIC_Main_Polarisation_H
#define AMEGIC_Main_Polarisation_H

#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"

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
    double Spin_Average(int,ATOOLS::Flavour*);
    int    Massless_Vectors(int,ATOOLS::Flavour*);
    int    Massive_Vectors(int,ATOOLS::Flavour*);
    void   Add_Extern_Polarisations(Basic_Sfuncs*,ATOOLS::Flavour*,Helicity*);
    void   Attach(int,ATOOLS::Flavour*);
    void   Reset_Gauge_Vectors(int,ATOOLS::Vec4D*,ATOOLS::Vec4D);
    void   Set_Gauge_Vectors(int,ATOOLS::Vec4D*,ATOOLS::Vec4D);
    double Massless_Norm(int,ATOOLS::Flavour*,Basic_Sfuncs*);
    double Massive_Norm() { return Mass_Norm; }
    void   Replace_Numbers(int,ATOOLS::Flavour*,Single_Amplitude*);
  };
}
#endif




