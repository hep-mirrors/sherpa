#ifndef Lund_Fortran_Interface_h
#define Lund_Fortran_Interface_h

#include "Blob_List.H"
#include "Particle_List.H"
#include "Particle.H"

namespace SHERPA {
  class Lund_Fortran_Interface {
    bool  pyinit;
    double a,b,sigma;

    static const int maxentries = 2000;
    double  * phep, * vhep;
    int     * jmohep, * jdahep, * isthep, * idhep;

    double  * pjet, * xjet;
    int       nk, * kfjet, * mothers, * daughters;

    void Init();
    void AddPartonToString(ATOOLS::Particle *,int &);
    void FillPrimaryHadronsInBlob(ATOOLS::Blob *,ATOOLS::Blob_List *,
				  ATOOLS::Particle_List * = NULL);
    void FillSecondaryHadronsInBlob(ATOOLS::Blob *,ATOOLS::Blob_List *,int,int,
				    ATOOLS::Particle_List * = NULL);
  public:
    Lund_Fortran_Interface(double,double,double);
    ~Lund_Fortran_Interface();
    bool Hadronize(ATOOLS::Blob *,ATOOLS::Blob_List *,
		   ATOOLS::Particle_List * = NULL);
  };

  extern Lund_Fortran_Interface LFInterface;
}

#endif
