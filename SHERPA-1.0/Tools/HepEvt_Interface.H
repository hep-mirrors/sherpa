#ifndef HepEvt_Interface_H
#define HepEvt_Interface_H

#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Flavour.H"
#include "Vector.H"

namespace SHERPA {
  class HepEvt_Interface {
    static const int maxentries = 4000;
    
    int m_evtnumber;

    std::map<ATOOLS::Particle *,int> m_connect;
    double * phep,   * vhep;
    int    * jmohep, * jdahep,  * isthep, * idhep;
    int      nhep;
    
    void ISBlobs2HepEvt(ATOOLS::Blob_List *,int &);
    void HardBlob2HepEvt(ATOOLS::Blob_List *,int &);
    void FSBlobs2HepEvt(ATOOLS::Blob_List *,int &);
    void FragmentationBlob2HepEvt(ATOOLS::Blob_List *,int &);
    void HadronDecayBlobs2HepEvt(ATOOLS::Blob_List *,int &);

    void EstablishRelations(ATOOLS::Particle *,ATOOLS::Blob *);  
    void EstablishRelations(ATOOLS::Blob *);  
    void Particle2HepEvt(ATOOLS::Particle *,int &);
    void String2HepEvt(ATOOLS::Blob *,int &);
  public:
    HepEvt_Interface();
    ~HepEvt_Interface();
    
    void Sherpa2HepEvt(ATOOLS::Blob_List *);

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

