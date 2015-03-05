#ifndef PHOTONS_Main_Photons_H
#define PHOTONS_Main_Photons_H

#include "MODEL/Main/Running_AlphaQED.H"
#include "PHOTONS++/Main/Dipole_Type.H"

namespace ATOOLS {
  class Blob;
  class Data_Reader;
#ifdef PHOTONS_DEBUG
  class Histogram_2D;
#endif
}

namespace PHOTONS {
  class Photons {
  private:
    std::string                 m_name;
    bool                        m_photonsadded;
    bool                        m_success;
  public:
    Photons(ATOOLS::Data_Reader* reader=NULL);
    Photons();
    ~Photons();
    bool AddRadiation(ATOOLS::Blob *);
    bool CheckStateBeforeTreatment(ATOOLS::Blob *);

    inline bool        DoneSuccessfully() { return m_success; }
    inline bool        AddedAnything()    { return m_photonsadded; }
    inline std::string Name()             { return m_name; }

    static double               s_alpha;
    static bool                 s_userunningparameters;
    static inline void ResetAlphaQED()
    { s_alpha=MODEL::aqed->AqedThomson(); }
    static inline void SetAlphaQED(const double& scale)
    { if (s_userunningparameters) s_alpha=(*MODEL::aqed)(scale); }

    static int                  s_mode;
    static bool                 s_useme;
    static double               s_ircutoff;
    static double               s_uvcutoff;
    static int                  s_ircutoffframe;
    static double               s_accu;
    static int                  s_nmax;
    static int                  s_nmin;
    static double               s_drcut;
    static bool                 s_strict;
    static double               s_increasemaxweight;
    static double               s_reducemaxenergy;
    static bool                 s_checkfirst;
    static int                  s_ffrecscheme;
    static int                  s_firecscheme;
#ifdef PHOTONS_DEBUG
    static std::string          s_histo_base_name;
    static ATOOLS::Histogram_2D s_histo_dipole;
    static ATOOLS::Histogram_2D s_histo_jacobianM;
    static ATOOLS::Histogram_2D s_histo_jacobianL;
    static ATOOLS::Histogram_2D s_histo_higher;
    static ATOOLS::Histogram_2D s_histo_yfs;
    static ATOOLS::Histogram_2D s_histo_total;
    static ATOOLS::Histogram_2D s_histo_t_dipole;
    static ATOOLS::Histogram_2D s_histo_t_jacobianM;
    static ATOOLS::Histogram_2D s_histo_t_jacobianL;
    static ATOOLS::Histogram_2D s_histo_t_higher;
    static ATOOLS::Histogram_2D s_histo_t_yfs;
    static ATOOLS::Histogram_2D s_histo_t_total;
#endif
  };


  

  /*!
    \namespace PHOTONS
    \brief namespace of the PHOTONS++ subprogramm
  */

  /*!
    \file Photons.H
    \brief Contains the class PHOTONS::Photons.
  */

  /*! 
    \class Photons
    \brief The Main Stearing Class for YFS Radiation in Sherpa
  */
  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Photons
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var bool Photons::m_analyse
    \brief bool given to constructor to switch treatment on or off
  */

  /*!
    \var bool Photons::m_success
    \brief value filled by Define_Dipole class, specifies whether treatment was successful
  */

  /*!
    \var bool Photons::m_photonsadded
    \brief value filled by Define_Dipole class, tells whether treatment added any photons
  */

  /*!
    \var static int Photons::s_mode
    \brief contains the mode of operation of Photons set in Fragmentation.dat
  */

  /*!
    \var static bool Photons::s_useme
    \brief if true, ME corrections will be enabled, set in Fragmentation.dat
  */ 

  /*!
    \var static double Photons::s_ircutoff
    \brief contains the infrared cut-off, set in Fragmentation.dat
  */
  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Photons
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Photons::Photons(ATOOLS::Data_Read * reader=NULL,bool=false)
    \brief sets the all variables to their values defined in Fragmentation.dat
  */

  /*!
    \fn Photons::Photons(bool=false)
    \brief initialises all variables to their default values
  */

  /*!
    \fn bool Photons::AddRadiation(Blob *)
    \brief hands blob to Define_Dipole class if it is supposed to be treated, returns m_photonsadded
  */

  /*!
    \fn bool Photons::DoneSuccessfully()
    \brief returns m_success
  */
} // end namespace PHOTONS
#endif
