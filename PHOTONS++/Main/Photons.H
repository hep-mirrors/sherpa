#ifndef PHOTONS_Main_Photons_H
#define PHOTONS_Main_Photons_H
#include "ATOOLS/Org/CXXFLAGS.H"


#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"
#include "PHOTONS++/Main/Define_Dipole.H"

using namespace ATOOLS;

namespace ATOOLS { class Data_Reader; }

namespace PHOTONS {
  class Photons {
  private:
    bool                        m_analyse;
    bool                        m_photonsadded;
    bool                        m_success;
  public:
    Photons(ATOOLS::Data_Reader* reader=NULL,bool=false);
    Photons(bool=false);
    ~Photons();
    bool AddRadiation(Blob *);
    bool DoneSuccessfully() { return m_success; }

    static int                  s_mode;
    static bool                 s_useme;
    static double               s_ircutoff;
    static double               s_accu;
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
