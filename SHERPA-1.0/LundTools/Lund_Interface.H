#ifndef Lund_Interface_H
#define Lund_Interface_H

#include "HepEvt_Interface.H"
#include "Particle_List.H"
#include "Blob_List.H"
#include "ISR_Handler.H"
#include <list>
#include <map>

namespace SHERPA {

  class Lund_Interface {
  private:
    static size_t s_errors, s_maxerrors;

    std::string m_path, m_file, m_outfile;
    size_t m_maxtrials;

    ATOOLS::HepEvt_Interface *p_hepevt;
    static PDF::ISR_Handler  *s_isrhandler;

    long int m_events, m_curfile, m_evtsperfile;
    bool     m_compress, m_writeout;

    std::list<ATOOLS::Particle_List *> m_partlists;
    double *p_phep, *p_vhep;
    int    *p_jmohep, *p_jdahep;

    static ATOOLS::Blob_List *s_bloblist;

    static bool s_exportas, s_exportpdf;

    std::map<ATOOLS::Particle *,int> m_secondarymap;
	
    void NextFile(const bool newfile);

    // preparing and doing the fragmentation
    void Reset();
    bool ExtractSinglets(ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *pl);
    int  PrepareFragmentationBlob(ATOOLS::Blob *blob);
    bool StringFragmentation(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
			     ATOOLS::Particle_List *pl,int nhep); 
    // tools for extracting hadrons
    void AddPartonToString(ATOOLS::Particle *particle,int &nhep);
    void FillPrimaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
				  ATOOLS::Particle_List *pl=NULL);
    void FillPrimaryTauInBlob(int,ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
			      ATOOLS::Particle_List *pl=NULL);
    void FillSecondaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
				    int daughter1,int daughter2,
				    ATOOLS::Particle_List *pl=NULL);
  public:

    Lund_Interface(std::string _m_path,std::string _m_file,bool sherpa);
    ~Lund_Interface();

    bool        Hadronize(ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *pl=NULL);
    bool        FillDecay(ATOOLS::Particle * part,
			  ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *pl=NULL);
    bool        FindDecay( ATOOLS::Particle * part );
    void        PerformAllDecays(ATOOLS::Blob *);
    void        PerformDecay( ATOOLS::Particle * part, 
			      ATOOLS::Blob_List * blob_list, 
			      ATOOLS::Particle_List * part_list );
    static void Error(const int error);
    bool        OneEvent(ATOOLS::Blob_List *const bloblist,double &weight);

    // inline functions
    inline static void SetISRHandler(PDF::ISR_Handler *const handler) { s_isrhandler=handler; }
    inline static PDF::ISR_Handler *ISRHandler()                      { return s_isrhandler; }

    inline static bool ExportAS()  { return s_exportas;  }
    inline static bool ExportPDF() { return s_exportpdf; }

    void SwitchOfDecays(ATOOLS::kf::code);
    void LeaveOpenOnlyDC( int idc ); 
  };// end of class Lund_Interface

}// end of namespace SHERPA

#endif
