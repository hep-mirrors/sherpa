#ifndef PHOTONS_PhaseSpace_Generate_Multipole_Photon_Angle_H
#define PHOTONS_PhaseSpace_Generate_Multipole_Photon_Angle_H

#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {
  class Particle;
  typedef std::vector<Particle*> Particle_Vector;
}

namespace PHOTONS {
  struct IdPairNbar;
  typedef std::vector<IdPairNbar> IdPairNbarVector;

  class Generate_Multipole_Photon_Angle {
    private:
      ATOOLS::Particle_Vector m_dipole;
      double              m_posnbar;
      IdPairNbarVector    m_nbars, m_posnbars;
      double              m_theta;
      double              m_phi;

      double  CalculateBeta(const ATOOLS::Vec4D&);
      void    GenerateMultipoleAngle();
      double  CalculateWeightByThetaPhi(const double&, const double&);
      double  CalculateWeightByVector(const ATOOLS::Vec4D&);
      double  TiTj(const size_t&, const size_t&);

    public:
      Generate_Multipole_Photon_Angle(const ATOOLS::Particle_Vector&,
                                      const IdPairNbarVector&);
      ~Generate_Multipole_Photon_Angle();

      inline double  GetTheta()      { return m_theta; }
      inline double  GetPhi()        { return m_phi; }
  };


  

  /*!
    \file Generate_Multipole_Photon_Angle.H
    \brief contains the class Generate_Multipole_Photon_Angle
  */

  /*!
    \class Generate_Multipole_Photon_Angle
    \brief generates the photon angular distribution in multipoles
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Generate_Multipole_Photon_Angle
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /*!
    \var Particle_Vector Generate_Multipole_Photon_Angle::m_dipole
    \brief contains the multipole
  */

  /*!
    \var std::vector<double> Generate_Multipole_Photon_Angle::m_nbars
    \brief avarage photon number of all constituent dipoles of the multipole
  */

  /*!
    \var double Generate_Multipole_Photon_Angle::m_theta
    \brief photon angle \f$ \theta \f$
  */

  /*!
    \var doubel Generate_Multipole_Photon_Angle::m_phi
    \brief photon azimutal angle \f$ \varphi \f$
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Generate_Multipole_Photon_Angle
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Generate_Multipole_Photon_Angle::Generate_Multipole_Photon_Angle(Particle_Vector, std::vector<double>)
    \brief initialises and executes the calculations

    The constructor <tt>Generate_Multipole_Photon_Angle::Generate_Multipole_Photon_Angle(Particle_Vector dip, std::vector<double> nbars)</tt> has to receive the multipole and the \f$ \bar{n} \f$ 's of every individual dipole. The numbering convention hereby is again as in \f$ \sum\limits_{i<j} \f$ , as is listed in the switch in <tt> Generate_Multipole_Photon_Angle::IndexLookup(unsigned int)</tt>. Therefore, here m_nbar \f$ =\sum|\bar{n}_k| \f$ and not \f$ \bar{n}=\sum\bar{n}_k \f$ as everywhere else.
  */

  /*!
    \fn double Generate_Multipole_Photon_Angle::CalculateBeta(Vec4D)
    \brief calculates \f$ \beta \f$ of a 4-vector
  */

  /*!
    \fn void Generate_Multipole_Photon_Angle::GenerateMultipoleAngle()
    \brief generates \f$ \theta \f$ and \f$ \varphi \f$

    One dipole is generated with probability \f$ \frac{|\bar{n}_k|}{\sum|\bar{n}_k|} \f$ and then the class Generate_Dipole_Photon_Angle is invoked to generate the photon angles
  */

  /*!
    \fn double Generate_Multipole_Photon_Angle::GetTheta()
    \brief returns m_theta ( \f$ \theta \f$ )
  */

  /*!
    \fn double Generate_Multipole_Photon_Angle::GetPhi()
    \brief returns m_phi ( \f$ \varphi \f$ )
  */

}

#endif
