#define IO_HepEvt_H
#ifdef  IO_HepEvt_H

#include "Parton_List.H"
#include "Blob_List.H"

namespace SHERPA {
  class IO_HepEvt {
    static const int maxentries = 2000;
    double * phep, * vhep;
    int    * jmohep, * jdahep, * isthep, * idhep;
    int    nhep;

    APHYTOOLS::Parton_List mypl;

    void MakeColorChain();
  public:
    IO_HepEvt();
    ~IO_HepEvt();
    void Blob2HepEvt(int, int mode =0);
    void Blob2HepEvt(APHYTOOLS::Blob *, int mode =0);
    void Blob2HepEvt(APHYTOOLS::Blob_List *, int mode =0);
    void Parton2HepEvt(APHYTOOLS::Parton *, int mode =0);
    void Parton2HepEvt(APHYTOOLS::Parton_List *, int mode =0);   
  };

}

#endif

