#ifndef IO_HepEvt_H
#define IO_HepEvt_H

#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Flavour.H"
#include "Vector.H"

namespace SHERPA {
  class IO_HepEvt {
    static const int maxentries = 2000;
    
    double * phep,   * vhep;
    int    * jmohep, * jdahep, * isthep, * idhep;
    int      nhep;
    
    void Parton2HepEvt(APHYTOOLS::Parton *,int *,int *,int);
    void MakeColorChain();
  public:
    IO_HepEvt();
    ~IO_HepEvt();
    
    void Blobs2HepEvt(APHYTOOLS::Blob_List *, std::string, int &);
    void Blob2HepEvt(APHYTOOLS::Blob *, int);
    void Entry2HepEvt(APHYTOOLS::Flavour,AMATOOLS::Vec4D,AMATOOLS::Vec4D,int);
    void Parton2HepEvt(APHYTOOLS::Parton *,int);

    int      Nhep()   { return nhep; }
    double * Phep()   { return phep; }
    double * Vhep()   { return vhep; }
    int    * Jmohep() { return jmohep; }
    int    * Jdahep() { return jdahep; }
    int    * Isthep() { return isthep; }
    int    * Idhep()  { return idhep; }
  };
}

#endif

