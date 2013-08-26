#ifndef PHOTONS_PhaseSpace_Avarage_Photon_Number_H
#define PHOTONS_PhaseSpace_Avarage_Photon_Number_H

#include "ATOOLS/Math/Vector.H"
#include "PHOTONS++/Main/Dipole_Type.H"

namespace PHOTONS {
  struct IdPairNbar;
  typedef std::vector<IdPairNbar> IdPairNbarVector;

  class Avarage_Photon_Number {
    private:
      double              m_omegaMax;
      double              m_omegaMin;
      ATOOLS::Particle_Vector m_dipole;

      double              m_nbar;
      IdPairNbarVector    m_nbars;

      void    CalculateAvaragePhotonNumber();
      double  CalculateBeta(const ATOOLS::Vec4D&);
      double  InterferenceTerm(const double&, const double&,
                               const double&, const double&);

    public:
      Avarage_Photon_Number(const ATOOLS::Particle_Vector&,
                            const double&, const double&);
      ~Avarage_Photon_Number();

      const double            GetNBar()  const { return m_nbar; }
      const IdPairNbarVector& GetNBars() const { return m_nbars; }

  };


  

  /*!
    \file Avarage_Photon_Number.H
    \brief contains the class Avarage_Photon_Number
  */

  /*!
    \class Avarage_Photon_Number
    \brief calculates \f$\bar{n}\f$ of the multipole
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member variables for Avarage_Photon_Number
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Avarage_Photon_Number::m_omegaMax
    \brief maximum photon energy allowed kinematically
  */

  /*!
    \var double Avarage_Photon_Number::m_omegaMin
    \brief minimum photon energy (infrared cut-off)
  */

  /*!
    \var int Avarage_Photon_Number::m_number1
    \brief number of sub-divisions for integration region one
  */

  /*!
    \var int Avarage_Photon_Number::m_number2
    \brief number of sub-divisions for integration region two
  */

  /*!
    \var int Avarage_Photon_Number::m_number3
    \brief number of sub-divisions for integration region three
  */

  /*!
    \var Particle_Vector Avarage_Photon_Number::m_dipole
    \brief contains all charged particles of the blob forming the multipole
  */

  /*!
    \var double Avarage_Photon_Number::m_nbar
    \brief the avarage photon number of the multipole
  */

  /*!
    \var std::vector<double> Avarage_Photon_Number::m_nbars
    \brief the avarage photon numbers of all individual constituent dipoles
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member methods for Avarage_Photon_Number
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Avarage_Photon_Number::Avarage_Photon_Number(Particle_Vector, double, double)
    \brief initialises and executes the calculations

    The constructor <tt>Avarage_Photon_Number::Avarage_Photon_Number(Particle_Vector dip, double wmax, double wmin)</tt> has to receive a Particle_Vector 
    containing all charged particles and the upper and lower 
    limit of the photon energy. The numbers of sub-divisions 
    for the different regions of the \f$\varphi\f$-integration 
    are set to
    - m_number1 = 100
    - m_number2 = 10000
    - m_number3 = 100
    .
    They are optimised for highly relativistic massive particles, 
    thus being somewhat over the top for very massive or less 
    rapid particles.
    Finally, <tt>Avarage_Photon_Number::CalculateAvaragePhotonNumber()</tt> 
    is called to calculate \f$\bar{n}\f$.
  */

  /*!
    \fn void Avarage_Photon_Number::CalculateAvaragePhotonNumber()
    \brief stears and combines the calculation of the different terms for \f$\bar{n}\f$

    \f$\bar{n}\f$ is calculated via
    \f[
      \bar{n} = \int\frac{d^3k}{k}\Theta(\Omega)\tilde{S}(k)
    \f]
    with
    \f[
      \tilde{S}(k) 
      = \sum_{i<j}\tilde{S}_{ij}(k) -\frac{\alpha}{4\pi^2}
          \sum_{i<j}Z_iZ_j\theta_i\theta_j
            \left(\frac{p_i}{(p_i\cdot k)}-\frac{p_j}{(p_j\cdot k)}\right)^2
    \f]
  */

  /*!
    \fn double Avarage_Photon_Number::CalculateBeta(Vec4D)
    \brief calculates \f$\beta\f$ of 4-vector
  */

  /*!
    \fn double Avarage_Photon_Number::Function(double, double, double, double)
    \brief calculates the value of the \f$\theta\f$-integrated unit sphere part of \f$\tilde{S}_{ij}\f$ at \f$\varphi\f$

    The Function is
    \f[
      missing
    \f]
    The passed values for <tt>Avarage_Photon_Number::Function(double bi, double bj, double alpha, double phi)</tt> are
    - \f$ \beta_i \f$
    - \f$ \beta_j \f$
    - \f$ \alpha = \angle(\vec{p}_i,\vec{p_j}) \f$ in the multipole-CMS
    - \f$ \varphi \f$
    .
  */

  /*!
    \fn double Avarage_Photon_Number::IntegralOne(double, double, double)
    \brief numerically integrates integration regions one and two (\f$\varphi\in[0,\frac{\pi}{2}]\f$)

    The passed values for <tt>Avarage_Photon_Number::IntegralOne(double bi, double bj, double alpha)</tt> are
    - \f$ \beta_i \f$
    - \f$ \beta_j \f$
    - \f$ \alpha = \angle(\vec{p}_i,\vec{p_j})\f$ in the multipole-CMS
    .
    Region one comprises of \f$ \varphi\in[0,a]\f$ and region two of 
    \f$ \varphi\in[a,\frac{\pi}{2}]\f$ with 
    \f$ a=1.5 - \frac{1}{50}\ln(\beta_i\beta_j) \f$ determined empirically.
  */

  /*!
    \fn double Avarage_Photon_Number::IntegralTwo(double, double, double)
    \brief numerically integrates integration region three (\f$\varphi\in[\pi,\frac{3\pi}{2}]\f$)

    The passed values for <tt>Avarage_Photon_Number::IntegralTwo(double bi, double bj, double alpha)</tt> are
    - \f$\beta_i\f$
    - \f$\beta_j\f$
    - \f$\alpha = \angle(\vec{p}_i,\vec{p_j})\f$ in the multipole-CMS
    .
  */

  /*!
    \fn double Avarage_Photon_Number::IntegrateInterferenceTerm(double, double, double)
    \brief integrates the term ~\f$\sin\alpha\f$, thus calls <tt>IntegralOne/Two(double,double,double)</tt>
  */

  /*!
    \fn double Avarage_Photon_Number::GetNBar()
    \brief returns m_nbar
  */

  /*!
    \fn std::vector<double> Avarage_Photon_Number::GetNBars()
    \brief returns m_nbars
  */

}

#endif
