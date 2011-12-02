#ifndef SHERPA_LundTools_Lund_Interface_H
#define SHERPA_LundTools_Lund_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Exception.H"
#include <list>
#include <map>

namespace ATOOLS {
  class Blob;
  class Blob_List;
  class Particle;
}

namespace SHERPA {
  class HepEvt_Interface;

  class Lund_Interface :
      public ATOOLS::Terminator_Object,
      public ATOOLS::Mass_Selector {
  private:
    static size_t s_errors, s_maxerrors;

    std::string m_path, m_file, m_outfile;
    size_t m_maxtrials;

    HepEvt_Interface *p_hepevt;

    long int m_events, m_curfile, m_evtsperfile;
    bool     m_compress, m_writeout;


    double *p_phep, *p_vhep;
    int    *p_jmohep, *p_jdahep;

    static ATOOLS::Blob_List *s_bloblist;

    static bool s_exportas;

    static int    * s_saved_mrpy;
    static double * s_saved_rrpy;

    void NextFile(const bool newfile);

    // preparing and doing the fragmentation
    ATOOLS::Return_Value::code ExtractSinglets(ATOOLS::Blob_List *);
    int                        PrepareFragmentationBlob(ATOOLS::Blob *);
    ATOOLS::Return_Value::code StringFragmentation(ATOOLS::Blob *,ATOOLS::Blob_List *,int); 
    // tools for extracting hadrons
    void AddPartonToString(ATOOLS::Particle *,int &);
    void FillFragmentationBlob(ATOOLS::Blob *);
    void FillOutgoingParticlesInBlob(ATOOLS::Blob *);

    void RestoreStatus();
    bool ReadInStatus(const std::string &path);
    void ReadInStatus(const std::string &file, int mode);
    void WriteOutStatus(const std::string &file);
    void PrepareTerminate();
  public:

    Lund_Interface(std::string _m_path,std::string _m_file);
    ~Lund_Interface();

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob_List *);
    ATOOLS::Return_Value::code PerformDecay(ATOOLS::Blob *);
    static void Error(const int error);

    ATOOLS::Return_Value::code OneEvent(ATOOLS::Blob_List *const,double &);

    bool IsAllowedDecay(kf_code);

    inline static bool ExportAS()  { return s_exportas;  }

    void SwitchOffDecays(kf_code kfc);
    void AdjustProperties(ATOOLS::Flavour flav);
    void SwitchOffMassSmearing();
    double GenerateMass(ATOOLS::Flavour flav, double min, double max);
    static void SaveStatus();
    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }

    ATOOLS::Flavour IdhepToSherpa(long int idhep);
    long int SherpaToIdhep(const ATOOLS::Flavour& flav);
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
    \fn Lund_Interface::PerformDecay( ATOOLS::Blob * )
    \brief Performs the decay of the incoming particle of <var>blob</var>.

    This is the main method for the decay of <var>blob</var>. Before decaying
    the blob, it checks a few conditions:
    -# blob_status::needs_hadrondecays (otherwise Return_Value::Error)
    -# 1 incoming, active particle (otherwise Return_Value::Error)
    -# mass of the incoming particle: peak-w_cut \< m* \< peak+w_cut where peak
       and w_cut are determined by \c PYDAT2/PMAS(KC,1) and \c PYDAT2/PMAS(KC,3)
       (otherwise Return_Value::Retry_Method)
    .

    If all these conditions are met, the decay is performed:
    -# the in-particle is set to unstable in the Pythia common block
    -# Pythia's \c PYDECY() routine is called for this decayer
    -# the resulting event is filled back into the blob
    -# the in-particle is set back to stable in Pythia
    .
    */
  /*!
    \fn Lund_Interface::SwitchOffDecays(kf_code kfc)
    \brief Sets the particle kfc stable

    This method tells Pythia to set the particle <var>kfc</var> stable. Might
    be obsolete, because the hadron decay handler will call the method
    SetAllStable.
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
