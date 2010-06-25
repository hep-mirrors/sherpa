#ifndef PHOTONS_Tools_YFS_Form_Factor_H
#define PHOTONS_Tools_YFS_Form_Factor_H

#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {
  class Particle;
  typedef std::vector<Particle* > Particle_Vector;
  class Gauss_Integrator;
  class Function_Base;
}

namespace PHOTONS {
  class YFS_Form_Factor {
    private:
      ATOOLS::Vec4D m_p1;
      ATOOLS::Vec4D m_p2;
      double      m_ks;
      double      m_m1;
      double      m_m2;
      double      m_Z1;
      double      m_Z2;
      double      m_x1;
      double      m_x2;
      double      m_xx1;
      double      m_xx2;
      double      m_Y;
      double      m_t1t2;

      ATOOLS::Function_Base *p_ig1, *p_ig2;
      ATOOLS::Gauss_Integrator *p_gi1, *p_gi2;

      double Y();
      double IntP1();
      double IntE();
      double IntP2();
      double GFunc(double);
      double IntG();
      double CalculateBeta(const ATOOLS::Vec4D&);

    public:
      YFS_Form_Factor(const ATOOLS::Particle_Vector&, const double&);
      YFS_Form_Factor(const ATOOLS::Particle*, const ATOOLS::Particle*,
                      const double&);
      ~YFS_Form_Factor();

      double G(double);

      inline const double Get() const { return m_Y; }

      inline const ATOOLS::Vec4D &P1() const  { return m_p1; }
      inline const ATOOLS::Vec4D &P2() const  { return m_p2; }
  };


  

  /*!
    \file YFS_Form_Factor.H
    \brief contains the class YFS_Form_Factor
  */

  /*!
    \class YFS_Form_Factor
    \brief calculates the YFS form factor of a multipole configuration given
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variable for YFS_Form_Factor
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Particle * YFS_Form_Factor::p_part1
    \brief particle 1 of the dipole
  */

  /*!
    \var Particle * YFS_Form_Factor::p_part2
    \brief particle 2 of the dipole
  */

  /*!
    \var Vec4D YFS_Form_Factor::m_p1
    \brief momentum of particle 1
  */

  /*!
    \var Vec4D YFS_Form_Factor::m_p2
    \brief momentum of particle 2
  */

  /*!
    \var double YFS_Form_Factor::m_ks
    \brief infrared cut-off
  */

  /*!
    \var double YFS_Form_Factor::m_m1
    \brief mass of particle 1
  */

  /*!
    \var double YFS_Form_Factor::m_m2
    \brief mass of particle 2
  */

  /*!
    \var double YFS_Form_Factor::m_Z1
    \brief charge of particle 1
  */

  /*!
    \var double YFS_Form_Factor::m_Z2
    \brief charge of particle 2
  */

  /*!
    \var double YFS_Form_Factor::m_x1
    \brief one of the roots of \f$ p_x^2 \f$
  */

  /*!
    \var double YFS_Form_Factor::m_x2
    \brief other root of \f$ p_x^2 \f$
  */

  /*!
    \var double YFS_Form_Factor::m_xx1
    \brief one of the roots of \f$ p_x^{'2} \f$
  */

  /*!
    \var double YFS_Form_Factor::m_xx2
    \brief other root of \f$ p_x^{'2} \f$
  */

  /*!
    \var double YFS_Form_Factor::m_Y
    \brief YFS form factor
  */

  /*!
    \var short int YFS_Form_Factor::m_t1t2
    \brief product of \f$ \theta_1\theta_2 \f$
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for YFS_Form_Factor
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn YFS_Form_Factor::YFS_Form_Factor(Particle_Vector, double)
    \brief multipole constructor

    The arguments to be passed are the multipole and the infrared cut-off (taken to be an isotropic energy cut-off in the same frame as the momenta of the multipole given). Determines the YFS form factor as
    \f[
      Y(\omega) = \sum_{i<j}Y_{ij}(\omega)
    \f]
    For the calculation of the individual dipole contributions the dipole constructor is used.
  */

  /*!
    \fn YFS_Form_Factor::YFS_Form_Factor(Particle *, Particle *, double)
    \brief dipole constructor, implements all variables and initiates the calculations

    The arguments to be passed are both dipole particles and the infrared cutoff (taken to be an isotropic energy cut-off in the same frame as the momenta of the multipole given). 

    The Particles are ordered such that \f$ E_2 \ge E_1 \f$, possible by the invariance of the form factor under the swapping of the labels of the particles.
  */

  /*!
    \fn double YFS_Form_Factor::Y()
    \brief calculates the YFS form factor

    \f[
      Y_{ij}(\omega) = 2\alpha\left(\mathcal{R}e B_{ij}+\tilde{B}_{ij}(\omega)\right)
    \f]
    \f[
      = -\frac{\alpha}{\pi}Z_1Z_2\theta_1\theta_2\left(\ln\frac{E_1E_2}{\omega^2}+\frac{1}{2}(p_1\cdot p_2)I_{p,1} - \frac{1}{2}(p_1\cdot p_2)I_{E} + \frac{1}{4}I_{p,2} + G(1) + G(-1) - (p_1 \cdot p_2)I_{G}\right)
    \f]
  */

  /*!
    \fn double YFS_Form_Factor::IntP1()
    \brief calculates the integral \f$ I_{p,1} \f$ (analytically)

    with \f$ I_{p,1} = \theta_1\theta_2\int\limits_{-1}^1dx\frac{\ln\frac{p_x^{'2}}{\lambda^2}}{p_x^{'2}} - \int\limits_{-1}^1dx\frac{\ln\frac{p_x^2}{\lambda^2}}{p_x^2} \f$
  */

  /*!
    \fn double YFS_Form_Factor::IntE()
    \brief calculates the integral \f$ I_E \f$ (analytically)

    with \f$ I_E = \int\limits_{-1}^1dx\frac{\ln\frac{E_x^2}{\omega^2}}{p_x^2} \f$
  */

  /*!
    \fn double YFS_Form_Factor::IntP2()
    \brief calculates the integral \f$ I_{p,2} \f$ (analytically)

    with \f$ I_{p,2} = \int\limits_{-1}^1dx\ln\frac{p_x^2}{m_1m_2} \f$
  */

  /*!
    \fn double YFS_Form_Factor::G(double)
    \brief calculates the function \f$ \tilde{G}(x) \f$ at the given value

    The value to be passed is the value at which \f$ \tilde{G}(x) = \frac{1-\beta_x}{2\beta_x}\ln\frac{1+\beta_x}{1-\beta_x}+\ln\frac{1+\beta_x}{2} \f$ shall be evaluated. The condition \f$ x \in [-1,1] \f$ must hold.
  */

  /*!
    \fn double YFS_Form_Factor::IntG()
    \brief calculates the integral \f$ I_G \f$ (numerically)

    with \f$ I_G = \int\limits_{-1}^1dx\frac{\tilde{G}(x)}{p_x^2} \f$. This integral cannot be solved analytically in a general frame.
  */

  /*!
    \fn double YFS_Form_Factor::CalculateBeta(Vec4D)
    \brief calculates \f$ \beta \f$ of a given 4-momentum
  */

  /*!
    \fn double YFS_Form_Factor::Get()
    \brief returns the YFS form factor
  */
}

#endif
