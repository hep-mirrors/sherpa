#ifndef PHOTONS_Tools_Generate_One_Photon_H
#define PHOTONS_Tools_Generate_One_Photon_H

#include "PHOTONS++/Main/Dipole_Type.H"

namespace PHOTONS {
  class Generate_One_Photon {
    private:
      double            m_beta1;
      double            m_beta2;
      double            m_theta;
      double            m_phi;
      double            m_w;
      double            m_delta;
      Dipole_Type::code m_dtype;
      ATOOLS::Particle *p_photon;

      void GeneratePhotonAngleMassless();
      void GeneratePhoton();

    public:
      Generate_One_Photon(const double&, const double&,
                          const double&, Dipole_Type::code);
      Generate_One_Photon(const ATOOLS::Particle_Vector&,
                          const std::vector<double>&,double, Dipole_Type::code);
      ~Generate_One_Photon();

      inline ATOOLS::Particle * GetPhoton() { return p_photon; }

      void SetMinimalPhotonAngle(const double&);
  };


  

  /*!
    \file Generate_One_Photon.H
    \brief contains the class Generate_One_Photon
  */

  /*!
    \class Generate_One_Photon
    \brief generates one photon according to the angular distribution of the multipole
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Generate_One_Photon
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Generate_One_Photon::m_beta1
    \brief \f$ \beta \f$ of particle 1 in a different mass dipole
  */

  /*!
    \var double Generate_One_Photon::m_beta2
    \brief \f$ \beta \f$ of particle 2 in a different mass dipole
  */

  /*!
    \var double Generate_One_Photon::m_theta
    \brief diced inclinational angle of the photon in a dipole
  */

  /*!
    \var double Generate_One_Photon::m_phi
    \brief diced azimutal angle of the photon in a dipole
  */

  /*!
    \var double Generate_One_Photon::m_w
    \brief energy for the photon
  */

  /*!
    \var double Generate_One_Photon::m_delta
    \brief collinear cut-off angle (needed for angular generation of massless particles)
  */

  /*!
    \var Dipole_Type::code Generate_One_Photon::m_dtype
    \brief specifies the Dipole_Type
  */

  /*!
    \var Particle * Generate_One_Photon::p_photon
    \brief the complete photon generated
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods fo Generate_One_Photon
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Generate_One_Photon::Generate_One_Photon(const double, const double, const double, Dipole_Type::code)
    \brief constructor for photon generation in dipoles

    The values passed are \f$ \beta_1,\beta_2 \f$, the energy of the photon and the Dipole_Type. For the dipole it is assumed that particle one is aligned along the \f$ +z \f$-axis for final-final dipoles and along the \f$ -z \f$-axis for initial-final dipoles, because in the this case, <tt>Dipole_FF</tt> fills the charged initial state into position 1. The constructor also executes the photon generation process.
  */

  /*!
    \fn Generate_One_Photon::Generate_One_Photon(Particle_Vector, std::vector<double>, double, Dipole_Type::code)
    \brief constructor for photon generation in multipoles

    The values passed are the multipole, all individual dipoles' \f$ \bar{n}_i \f$, the energy of the photon and the Dipole_Type. The constructor also executes the photon generation process wherein for the angular generation <tt>Generate_Multipole_Photon_Angle</tt> is called. It is also able to handle dipoles, but somewhat slower than the built-in generation.
  */

  /*!
    \fn void Generate_One_Photon::GeneratePhotonAngleMassless()
    \brief generates \f$ \theta \f$ in massless dipoles according to \f$ \tilde{S}_{ij} \f$
  */

  /*!
    \fn void Generate_One_Photon::GeneratePhoton()
    \brief builds a particle out of the generated angles and the energy
  */

  /*!
    \fn Particle * Generate_One_Photon::GetPhoton()
    \brief returns p_photon, the generated photon
  */

  /*!
    \fn void Generate_One_Photon::SetMinimalPhotonAngle(double)
    \brief changes the preset value of the collinear cut-off angle needed in the massless generation and regenerates the event
  */

}

#endif
