#ifndef HepEvt_Interface_H
#define HepEvt_Interface_H

#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Flavour.H"
#include "Vector.H"

namespace SHERPA {
  struct gtp {
    enum code {
      Unspecified  = 0,
      Sherpa       = 1,
      Herwig       = 2,
      Pythia       = 3
    };
  };

  class HepEvt_Interface {
    bool            m_io;
    int             m_mode;
    std::string     m_path, m_file;
    std::ofstream   m_outstream;
    std::ifstream * p_instream;

    ATOOLS::Particle_Int_Map m_convertS2H;
    ATOOLS::Int_Particle_Map m_convertH2S;

    int       m_evtnumber,m_nhep;
    int       m_filesize, m_evtcount;
    gtp::code m_generator;

    static const int s_maxentries = 4000;
    double * p_phep,   * p_vhep;
    int    * p_jmohep, * p_jdahep,  * p_isthep, * p_idhep;
    double   m_weight;
    
    void WriteFullHepEvt(int nevt);
    void WriteReducedHepEvt(int nevt);
    void ISBlobs2HepEvt(ATOOLS::Blob_List * const,int &);
    void HardBlob2HepEvt(ATOOLS::Blob_List * const,int &);
    void FSBlobs2HepEvt(ATOOLS::Blob_List * const,int &);
    void FragmentationBlob2HepEvt(ATOOLS::Blob_List * const,int &);
    void HadronDecayBlobs2HepEvt(ATOOLS::Blob_List * const,int &);

    void EstablishRelations(ATOOLS::Particle * const,ATOOLS::Blob * const);  
    void EstablishRelations(ATOOLS::Blob * const);  
    void Particle2HepEvt(ATOOLS::Particle * const,int &);
    void String2HepEvt(ATOOLS::Blob * const,int &);



    bool ReadHepEvt(ATOOLS::Blob_List * const);
    void HepEvt2Particle(const int,const int);
    bool ConstructBlobsFromHerwig(ATOOLS::Blob_List * const);
    bool ConstructBlobsFromPythia(ATOOLS::Blob_List * const);
    bool ConstructBlobs(ATOOLS::Blob_List * const);
    bool IdentifyBlobs(ATOOLS::Blob_List * const);
    bool DeleteObsolete();
  public:

    HepEvt_Interface(bool=true,int=11,
		     const std::string & =std::string(""),
		     const std::string & =std::string(""));
    ~HepEvt_Interface();
    
    int                  Nhep() const   { return m_nhep;   }
    const double * const Phep() const   { return p_phep;   }
    const double * const Vhep() const   { return p_vhep;   }
    const int    * const Jmohep() const { return p_jmohep; }
    const int    * const Jdahep() const { return p_jdahep; }
    const int    * const Isthep() const { return p_isthep; }
    const int    * const Idhep() const  { return p_idhep;  }
    double               Weight()       { return m_weight; }

    bool Sherpa2HepEvt(ATOOLS::Blob_List * const);
    bool HepEvt2Sherpa(ATOOLS::Blob_List * const);
  };
}

#endif

