#ifndef Lund_Fortran_Interface_h
#define Lund_Fortran_Interface_h

#include "Parton_List.H"
#include "Parton.H"

namespace SHERPA {
  class Lund_Fortran_Interface {
    bool  pyinit;
    double a,b,sigma;

    static const int maxentries = 2000;
    double * phep, * vhep;
    int    * jmohep, * jdahep, * isthep, * idhep;

    void Init();
    void AddPartonToString(APHYTOOLS::Parton *,int &);
    void FillHadronsInBlob(APHYTOOLS::Parton_List *,APHYTOOLS::Blob *);
  public:
    Lund_Fortran_Interface(double,double,double);
    ~Lund_Fortran_Interface();
    void Hadronize(APHYTOOLS::Parton_List *,APHYTOOLS::Blob *);
  };

  extern Lund_Fortran_Interface LFInterface;
}

#endif
