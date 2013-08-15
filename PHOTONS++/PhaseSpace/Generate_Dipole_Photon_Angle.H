#ifndef PHOTONS_PhaseSpace_Generate_Dipole_Photon_Angle_H
#define PHOTONS_PhaseSpace_Generate_Dipole_Photon_Angle_H

#include "ATOOLS/Math/Vector.H"

namespace PHOTONS {
  class Generate_Dipole_Photon_Angle {
    private:
      double              m_b1;
      double              m_b2;
      double              m_c;
      double              m_theta;
      double              m_phi;
      ATOOLS::Vec4D       m_dir;

      double  CalculateBeta(const ATOOLS::Vec4D&);
      void    GenerateDipoleAngle();
      void    GenerateNullVector();

    public:
      Generate_Dipole_Photon_Angle(ATOOLS::Vec4D, ATOOLS::Vec4D);
      Generate_Dipole_Photon_Angle(const double&, const double&);
      ~Generate_Dipole_Photon_Angle();

      inline double               GetCosTheta()   { return m_c; }
      inline double               GetTheta()      { return m_theta; }
      inline double               GetPhi()        { return m_phi; }
      inline const ATOOLS::Vec4D& GetVector()     { return m_dir; }
  };


  

  /*!
    \file Generate_Dipole_Photon_Angle.H
    \brief contains the class Generate_Dipole_Photon_Angle
  */

  /*!
    \class Generate_Dipole_Photon_Angle
    \brief generates the photon angular distribution in dipoles
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member variables for Generate_Dipole_Photon_Angle
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Generate_Dipole_Photon_Angle::m_b1
    \brief \f$ \beta_1 \f$
  */

  /*!
    \var double Generate_Dipole_Photon_Angle::m_b2
    \brief \f$ \beta_2 \f$
  */

  /*!
    \var double Generate_Dipole_Photon_Angle::m_c
    \brief \f$ c = \cos\theta \f$
  */

  /*!
    \var double Generate_Dipole_Photon_Angle::m_theta
    \brief \f$ \theta \f$
  */

  /*!
    \var double Generate_Dipole_Photon_Angle::m_phi
    \brief \f$ \varphi \f$
  */

  /*!
    \var Vec4D Generate_Dipole_Photon_Angle::m_dir
    \brief null vector of unit spatial length in direction \f$ (\theta,\varphi) \f$
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Generate_Dipole_Photon_Angle
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Generate_Dipole_Photon_Angle::Generate_Dipole_Photon_Angle(Vec4D, Vec4D)
    \brief generates dipole angles for two arbitrary timelike 4-vectors \f$ p_1 \f$ and \f$ p_2 \f$

    \f$ p_1 \f$ and \f$ p_2 \f$ are boosted in their CMS, there the photon angles are calculated and m_dir is generated. Finally, m_dir is boosted to the original system of \f$ p_1 \f$ and \f$ p_2 \f$ and \f$ \theta \f$ and \f$ \varphi \f$ are recalculated.

    This constructor is used by the Generate_Multipole_Photon_Angle class.
  */

  /*!
    \fn Generate_Dipole_Photon_Angle::Generate_Dipole_Photon_Angle(double, double)
    \brief generates dipole angles for two 4-vectors with \f$ \beta_1 \f$ and \f$ \beta_2 \f$ assumed to be in their CMS and aligned along the z-axis

    Both angles are calculated via <tt>GenerateDipoleAngle()</tt>. No null vector will be produced.

    This constructor is used by the Generate_One_Photon class.
  */

  /*!
    \fn double Generate_Dipole_Photon_Angle::CalculateBeta(Vec4D)
    \brief calculates \f$ \beta \f$ for a given 4-vector
  */

  /*!
    \fn void Generate_Dipole_Photon_Angle::GenerateDipoleAngle()
    \brief generates both photon angles

    Works in the dipole's CMS. \f$ \varphi \f$ is distributed uniformly, \f$ \theta \f$ according to the eikonal factor \f$ \tilde{S}_{ij} \f$ .
  */

  /*!
    \fn void Generate_Dipole_Photon_Angle::GenerateNullVector()
    \brief m_dir is generated

    This null vector can be Poincare transformed to any frame to have the photon angular configuration there. To get the full photon its energy/3-momentum simply has to be multiplied by the generated energy.
  */

  /*!
    \fn double Generate_Dipole_Photon_Angle::GetCosTheta()
    \brief returns m_c ( \f$ c = \cos\theta \f$ )
  */

  /*!
    \fn double Generate_Dipole_Photon_Angle::GetTheta()
    \brief returns m_theta ( \f$ \theta \f$ )
  */

  /*!
    \fn double Generate_Dipole_Photon_Angle::GetPhi()
    \brief returns m_phi ( \f$ \varphi \f$ )
  */

  /*!
    \fn Vec4D Generate_Dipole_Photon_Angle::GetVector()
    \brief returns m_dir
  */


}


// this class will take two four vectors and generate a null vector of unit 3D length which is distributed according to eikonal factor
// if two doubles (b1,b2) are given it assumed they are in their respective rest frame and then this vector is generated in that frame
#endif
