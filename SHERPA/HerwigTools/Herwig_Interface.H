#ifndef SHERPA_HerwigTools_Herwig_Interface_H
#define SHERPA_HerwigTools_Herwig_Interface_H

#include "SHERPA/Tools/HepEvt_Interface.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "PDF/Main/ISR_Handler.H"
#include <map>

namespace SHERPA {
  inline void MakeFortranString(char * output,std::string input,unsigned int length)
  {
    for (unsigned int i=0;i<length;++i)   output[i]=(char)32;
    for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
  }


  class Herwig_Interface {
  private:
    static size_t s_errors, s_maxerrors;

    std::string m_path, m_file, m_outfilename;

    ATOOLS::HepEvt_Interface *p_hepevt;
    static PDF::ISR_Handler  *s_isrhandler;

    long int m_events, m_curfile, m_evtsperfile;
    bool     m_compress, m_writeout;

    double * p_phep,   * p_vhep;
    int    * p_jmohep, * p_jdahep;
    static ATOOLS::Blob_List * s_bloblist;

    static bool s_exportas, s_exportpdf;

    void ReadInTheParameters();
    void NextFile(const bool newfile);
  public:
    Herwig_Interface(std::string _m_path,std::string _m_file,bool sherpa);
    ~Herwig_Interface();

    bool Initialize();
    static void Error(const int error);
    bool OneEvent(ATOOLS::Blob_List *const bloblist,double &weight);

   // inline functions
    inline static void SetISRHandler(PDF::ISR_Handler *const handler) { s_isrhandler=handler; }

    inline static PDF::ISR_Handler *ISRHandler() { return s_isrhandler; }

    inline static bool ExportAS()  { return s_exportas;  }
    inline static bool ExportPDF() { return s_exportpdf; }

    //    bool Treat(ATOOLS::Blob_List *, double &);
    //    void CleanUp();
  };

}// end of namespace SHERPA

#endif
