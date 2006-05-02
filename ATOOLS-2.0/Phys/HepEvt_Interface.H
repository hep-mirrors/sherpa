#ifndef HepEvt_Interface_H
#define HepEvt_Interface_H

#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Flavour.H"
#include "Vector.H"

namespace ATOOLS {
  struct gtp {
    enum code {
      Unspecified  = 0,
      Sherpa       = 1,
      Herwig       = 2,
      Pythia       = 3
    };
  };
  typedef std::map<int,std::pair<Particle *,bool> > Translation_Map;

  class HepEvt_Interface {
    bool            m_hadronized;
    bool            m_io;
    int             m_mode;
    std::string     m_path, m_file;
    std::ifstream * p_instream;
    std::ofstream * p_outstream;

    ATOOLS::Particle_Int_Map m_convertS2H;
    Translation_Map          m_convertH2S;

    int       m_evtnumber,m_nhep;
    long int  m_filesize, m_evtcount;
    gtp::code m_generator;

    static const int s_maxentries = 4000;
    double * p_phep,   * p_vhep;
    int    * p_jmohep, * p_jdahep,  * p_isthep, * p_idhep;
    double   m_weight, m_Q2, m_x1, m_x2;
    int      m_fl1, m_fl2;
    
    void WriteFullHepEvt(int nevt);
    void WriteFormatedHepEvt(int nevt);
    void WriteReducedHepEvt(int nevt);
    void WriteD0HepEvt(int nevt);
    void ISBlobs2HepEvt(ATOOLS::Blob_List * const,int &);
    void HardBlob2HepEvt(ATOOLS::Blob_List * const,int &);
    void FSBlobs2HepEvt(ATOOLS::Blob_List * const,int &);
    void FragmentationBlob2HepEvt(ATOOLS::Blob_List * const,int &);
    void HadronDecayBlobs2HepEvt(ATOOLS::Blob_List * const,int &);

    void EstablishRelations(ATOOLS::Particle * const,ATOOLS::Blob * const);  
    void EstablishRelations(ATOOLS::Blob * const);  
    void Particle2HepEvt(ATOOLS::Particle * const,int &);
    void String2HepEvt(ATOOLS::Blob * const,int &);
    
    
    void ChangeOutStream();

    void ReadHepEvt(ATOOLS::Blob_List * const);
    void OpenNewHepEvtFile();
    void HepEvt2Particle(const int);
    bool ConstructBlobsFromHerwig(ATOOLS::Blob_List * const);
    void FollowDaughters(const int & i,ATOOLS::Blob * hadron,ATOOLS::Blob * clusters);
    bool ConstructBlobsFromPythia(ATOOLS::Blob_List * const);
    bool ConstructBlobs(ATOOLS::Blob_List * const);
    bool IdentifyBlobs(ATOOLS::Blob_List * const);
    void DeleteObsolete(const int mode);
  public:
    HepEvt_Interface(int=0);
    HepEvt_Interface(gtp::code=gtp::Unspecified);
    HepEvt_Interface(bool=true,int=11,
		     const std::string & =std::string(""),
		     const std::string & =std::string(""),
		     const int _filesize=1000);
    ~HepEvt_Interface();
    
    int                  Nhep() const   { return m_nhep;   }
    const double * const Phep() const   { return p_phep;   }
    const double * const Vhep() const   { return p_vhep;   }
    const int    * const Jmohep() const { return p_jmohep; }
    const int    * const Jdahep() const { return p_jdahep; }
    const int    * const Isthep() const { return p_isthep; }
    const int    * const Idhep() const  { return p_idhep;  }
    double               Weight()       { return m_weight; }

    void SetNhep(int      _nhep)   { m_nhep = _nhep; }
    void SetPhep(double * _phep)   { p_phep = _phep; }
    void SetVhep(double * _vhep)   { p_vhep = _vhep; }
    void SetJmohep(int  * _jmohep) { p_jmohep = _jmohep; }
    void SetJdahep(int  * _jdahep) { p_jdahep = _jdahep; }
    void SetIsthep(int  * _isthep) { p_isthep = _isthep; }
    void SetIdhep(int   * _idhep)  { p_idhep = _idhep; }
    void SetWeight(double _weight) { m_weight = _weight; }
    void SetQ2(double _Q2) { m_Q2 = _Q2; }
    void Setx1(double _x1) { m_x1 = _x1; }
    void Setx2(double _x2) { m_x2 = _x2; }
    void SetFl1(int   _fl1) { m_fl1 = _fl1; }
    void SetFl2(int   _fl2) { m_fl2 = _fl2; }

    std::ofstream * GetOutStream() { return p_outstream; }
    void ChangeOutStream(std::string &,long int=-1);
    void PrintHepEvtEvent(const int=0);

    bool Sherpa2HepEvt(ATOOLS::Blob_List * const);
    bool HepEvt2Sherpa(ATOOLS::Blob_List * const);
    inline void SetHadronized(const bool hadronized) { m_hadronized=hadronized; }
    inline bool Hadronized() const                   { return m_hadronized;     }
  };
}

#endif

