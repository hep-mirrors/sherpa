#ifndef PHOTONS_Main_Define_Dipole_H
#define PHOTONS_Main_Define_Dipole_H

#include "PHOTONS++/Main/Dipole_Type.H"

namespace ATOOLS {
  class Blob;
}

namespace PHOTONS {
  class Define_Dipole {
    private:
      bool                    m_success;
      bool                    m_photonsadded;
      ATOOLS::Blob *          p_blob;
      ATOOLS::Particle_Vector m_chargedinparticles;
      ATOOLS::Particle_Vector m_neutralinparticles;
      ATOOLS::Particle_Vector m_chargedoutparticles;
      ATOOLS::Particle_Vector m_neutraloutparticles;
      ATOOLS::Particle_Vector m_softphotons;
      Particle_Vector_Vector  m_pvv;
      Dipole_Type::code       m_dtype;

    public:
      Define_Dipole(ATOOLS::Blob *);
      ~Define_Dipole();

      void AddRadiation();

      inline bool DoneSuccessfully()             { return m_success; }
      inline bool AddedAnything()                { return m_photonsadded; }
  };


  

  /*!
    \file Define_Dipole.H
    \brief contains the class Define_Dipole
  */

  /*!
    \class Define_Dipole
    \brief classifies and coordinates the treatment of the multipoles given by Photons class
  */
  ///////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Define_Dipole
  ///////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var bool Define_Dipole::m_success
    \brief value filled by Dipole_FF/Dipole_FI classes, specifies whether treatment was successful
  */

  /*!
    \var bool Define_Dipole::m_photonsadded
    \brief value filled by Dipole_FF/Dipole_FI classes, tells whether treatment added any photons
  */

  /*!
    \var Blob * Define_Dipole::p_blob
    \brief pointer to blob passed from Photons class
  */

  /*!
    \var Particle_Vector Define_Dipole::m_chargedinparticles
    \brief contains all charged initial state particles of p_blob
  */

  /*!
    \var Particle_Vector Define_Dipole::m_neutralinparticles
    \brief contains all neutral initial state particles of p_blob
  */

  /*!
    \var Particle_Vector Define_Dipole::m_chargedoutparticles
    \brief contains all charged final state particles of p_blob
  */

  /*!
    \var Particle_Vector Define_Dipole::m_neutraloutparticles
    \brief contains all neutral final state particles of p_blob
  */

  /*!
    \var Particle_Vector Define_Dipole::m_softphotons
    \brief contains all photons created by Dipole_FF/Dipole_FI classes
  */

  /*!
    \var Particle_Vector_Vector Define_Dipole::m_pvv
    \brief contains ensemble of all particles of p_blob in respective Particle_Vector's
  */

  /*!
    \var Dipole_Type::code Define_Dipole::m_dtype
    \brief contains the dipole type of p_blob
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods of Define_Dipole
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Define_Dipole::Define_Dipole(Blob *)
    \brief initialises the main parameters for the treatment of the blob

    In particular:
    - classifies all particles of the blob into the respective Particle_Vector
    - determines the Dipole_Type of the blob
    - creates m_pvv to hand over to Dipole_FF/Dipole_FI classes
    .
  */

  /*!
    \fn void Define_Dipole::AddRadiation()
    \brief gives m_pvv to the correct treatment class, reconstructs the blob after treatment

    In more detail:
    - hands m_pvv over to the Dipole_FF/Dipole_FI classes according to the 
      Dipole_Type determined in constructor
    - if treatment was successful, i.e. momentum conservation is fulfilled, 
      the p_blob will be reconstructed with the corrected particle momenta 
      and all photons
    - values of m_success and m_addedphotons are assigned
    .
  */

  /*!
    \fn bool Define_Dipole::DoneSuccessfully()
    \brief returns m_success
  */

  /*!
    \fn bool Define_Dipole::AddedAnything()
    \brief returns m_photonsadded
  */


}
#endif
