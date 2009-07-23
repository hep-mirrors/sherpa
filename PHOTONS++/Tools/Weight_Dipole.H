#ifndef PHOTONS_Tools_Weight_Dipole_H
#define PHOTONS_Tools_Weight_Dipole_H

#include "PHOTONS++/Tools/Weight_Base.H"


using namespace ATOOLS;

namespace PHOTONS {
  class Weight_Dipole: public Weight_Base {
    private:
      Particle_Vector     m_olddipole;
      Particle_Vector     m_newdipole;
      Particle_Vector     m_softphotons;

      virtual void CalculateWeight();
      virtual void CalculateMax();
      double SMod(const Vec4D&, const Vec4D&, const Vec4D&);

    public:
      Weight_Dipole(const Particle_Vector&, const Particle_Vector&,
                    const Particle_Vector&, Dipole_Type::code);
      virtual ~Weight_Dipole();
  };


  

  /*!
    \file Weight_Dipole.H
    \brief contains the class Weight_Dipole
  */

  /*!
    \class Weight_Dipole
    \brief calculates \f$ W_{dipole} \f$ of the event given
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the variables for Weight_Dipole
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Particle_Vector Weight_Dipole::m_olddipole
    \brief contains all charged particles before the treatment (\f$ K=0 \f$ event)
  */

  /*!
    \var Particle_Vector Weight_Dipole::m_newdipole
    \brief contains all charged particles after the treatment
  */

  /*!
    \var Particle_Vector Weight_Dipole::m_softphotons
    \brief contains all bremsstrahlung photons generated in the event
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member methods for Weight_Dipole
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Dipole::Weight_Dipole(const Particle_Vector, const Particle_Vector, const Particle_Vector, Dipole_Type::code)
    \brief initialises all variales and calculates the weights

    The values to be passed are the multipole before the treatment, after the treatment, the generated photons and the Dipole_Type.

    Calls <tt>Weight_Dipole::CalculateWeight()</tt> and <tt>Weight_Dipole::CalculateMax()</tt> to calculate the weight and its maximum. The latter value is the contribution to the maximum total weight as well as the maximum of \f$ W_{dipole} \f$ by itself.
  */

  /*!
    \fn void Weight_Dipole::CalculateWeight()
    \brief calculates the weight of the event

    The weight is \f$ W_{dipole} = \prod_{i=1}^{n_\gamma}\frac{\tilde{S}(p,P,k_i)}{\tilde{S}(p,Q,k_i)} \f$

    with \f$ \tilde{S}(k) = \sum\limits_{i<j}\tilde{S}_{ij}(k) = \sum\limits_{i<j}Z_iZ_j\theta_i\theta_j\left(\frac{p_1}{(p_1\cdot k_i)}-\frac{p_2}{(p_2\cdot k_i)}\right)^2 \f$

    Calls <tt>Weight_Dipole::SMod(const Vec4D, const Vec4D, const Vec4D)</tt> to calculate \f$ \tilde{S} \f$.
  */

  /*!
    \fn void Weight_Dipole::CalculateMax()
    \brief calculates the maximum weight of the configuration

    The maximum weight is \f$ W_{dipole,max} = 1 \f$.
  */

  /*!
    \fn double Weight_Dipole::SMod(const Vec4D, const Vec4D, const Vec4D)
    \brief calculates the eikonal factor in \f$ \tilde{S}_{ij} \f$

    The arguments are \f$ p_1 \f$, \f$ p_2 \f$ and \f$ k_i \f$.
    The value returned is \f$ \left(\frac{p_1}{(p_1\cdot k_i)}-\frac{p_2}{(p_2\cdot k_i)}\right)^2 \f$.
  */

}

#endif
