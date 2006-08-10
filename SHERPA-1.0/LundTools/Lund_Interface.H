#ifndef Lund_Interface_H
#define Lund_Interface_H

#include "HepEvt_Interface.H"
#include "Particle_List.H"
#include "Blob_List.H"
#include "ISR_Handler.H"
#include "Return_Value.H"
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


    double *p_phep, *p_vhep;
    int    *p_jmohep, *p_jdahep;

    static ATOOLS::Blob_List *s_bloblist;

    static bool s_exportas, s_exportpdf;

    void NextFile(const bool newfile);

    // preparing and doing the fragmentation
    ATOOLS::Return_Value::code ExtractSinglets(ATOOLS::Blob_List *);
    int                        PrepareFragmentationBlob(ATOOLS::Blob *);
    bool                       StringFragmentation(ATOOLS::Blob *,ATOOLS::Blob_List *,int); 
    // tools for extracting hadrons
    void AddPartonToString(ATOOLS::Particle *,int &);
    void FillFragmentationBlob(ATOOLS::Blob *);
    void FillOutgoingParticlesInBlob(ATOOLS::Blob *);
  public:

    Lund_Interface(std::string _m_path,std::string _m_file,bool sherpa);
    ~Lund_Interface();

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob_List *);
    ATOOLS::Return_Value::code PerformDecay(ATOOLS::Blob *);
    static void Error(const int error);

    ATOOLS::Return_Value::code OneEvent(ATOOLS::Blob_List *const,double &);

    bool IsAllowedDecay(ATOOLS::kf::code);
    // inline functions
    inline static void SetISRHandler(PDF::ISR_Handler *const handler) { s_isrhandler=handler; }
    inline static PDF::ISR_Handler *ISRHandler()                      { return s_isrhandler; }

    inline static bool ExportAS()  { return s_exportas;  }
    inline static bool ExportPDF() { return s_exportpdf; }

    void SwitchOfDecays(ATOOLS::kf::code kfc);
  };// end of class Lund_Interface

  /*!
    \file Lund_Interface.H
    \brief Declares the class SHERPA::Lund_Interface
    */
  /*!
    \class Lund_Interface
    \brief Tools for calling Pythia routines
    */
  /*!
    \fn Lund_Interface::PerformDecay( ATOOLS::Particle * part, ATOOLS::Blob_List * blob_list, ATOOLS::Particle_List * part_list )
    \brief Performs the decay of the particle <var>part</var>

    This is the main method for the decay of the particle <var>part</var>
    -# It choses a decay channel  and its kinemtaic by calling Pythia common blocks
    -# A new blob is created.
    -# Treats every daughter self-recursively.
    .
    */
  /*!
    \fn Lund_Interface::SwitchOfDecays(ATOOLS::kf::code kfc)
    \brief Sets the particle kfc stable

    This method tells Pythia to set the particle <var>kfc</var> stable.
    */
  /*!
    \fn Lund_Interface::FillPrimaryTauInBlob(int,ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist, ATOOLS::Particle_List *pl=NULL)
    \brief Reconstructs the tau-lepton and pushs it through the Blob

    This method became necessery for the implementation of tau decays into HADRONS.
    Basically, it pushs a (anti-)tau through the Fragmentation blob such that it can be treated
    by Hadron_Decays::Treat.
    */
  /*!
    \fn Lund_Interface::FillPrimaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist, ATOOLS::Particle_List *pl=NULL)
    \brief Method responsible for the hadronization and filling the Fragmentation Blob

    It's this methods task to deal with the hadronization and fill the Fragmentation Blob such that
    the outgoing hadrons can be treated by Hadron_Decays::Treat. If there are however incoming tau-leptons 
    they need to be "pushed through the blob" by calling Lund_Interface::FillPrimaryTauInBlob.
    */
  /*!
    \fn Lund_Interface::FindDecay( ATOOLS::Particle * part )
    \brief Returns true of Pythia can cope with the decay of the particle
    */


}// end of namespace SHERPA

#endif
