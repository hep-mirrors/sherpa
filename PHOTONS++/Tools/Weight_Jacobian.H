#ifndef PHOTONS_Tools_Weight_Jacobian_H
#define PHOTONS_Tools_Weight_Jacobian_H

#include "PHOTONS++/Tools/Weight_Base.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"

namespace PHOTONS {
  class Weight_Jacobian: public Weight_Base {
    protected:
      Vec4D   CalculateMomentumSum(Particle_Vector);

      virtual void CalculateWeight() = 0;
      virtual void CalculateMax() = 0;

    public:
      Weight_Jacobian();
      virtual ~Weight_Jacobian();
  };


  class Weight_Jacobian_Lorentz: public Weight_Jacobian {
    private:
      double    m_M;
      double    m_QC0;
      double    m_QN0;
      double    m_PC0;
      double    m_PN0;
      double    m_K0;
      double    m_mMQ;
      double    m_mMP;

      virtual void CalculateWeight();
      virtual void CalculateMax();

    public:
      Weight_Jacobian_Lorentz(Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, Dipole_Type::code);
      virtual ~Weight_Jacobian_Lorentz();
  };


  class Weight_Jacobian_Mapping: public Weight_Jacobian {
    private:
      Particle_Vector       m_newdipole;
      Particle_Vector       m_olddipole;
      Particle_Vector       m_newspectator;
      Particle_Vector       m_oldspectator;
      double                m_M;
      double                m_u;
      Vec3D                 m_K;

      virtual void CalculateWeight();
      virtual void CalculateMax();

    public:
      Weight_Jacobian_Mapping(Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, Vec4D, double, double, Dipole_Type::code);
      virtual ~Weight_Jacobian_Mapping();
  };


  

  /*!
    \file Weight_Jacobian.H
    \brief contains the classes Weight_Jacobian and its daughters Weight_Jacobian_Lorentz and Weight_Jacobian_Mapping
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /*!
    \class Weight_Jacobian
    \brief base class for both Jacobian weights classes Weight_Jacobian_Lorentz and Weight_Jacobian_Mapping
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Weight_Jacobian
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Jacobian::Weight_Jacobian()
    \brief contructor of purely virtual class
  */

  /*!
    \fn Vec4D Weight_Jacobian::CalculateMomentumSum(Particle_Vector)
    \brief calculates the sum of the momenta of the particles in the Particle_Vector given
  */

  /*!
    \fn void Weight_Jacobian::CalculateWeight()
    \brief purely virtual in Weight_Jacobian
  */

  /*!
    \fn void Weight_Jacobian::CalculateMax()
    \brief purely virtual in Weight_Jacobian
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /*!
    \class Weight_Jacobian_Lorentz
    \brief calculates \f$ W_{J,L} \f$, the weight due to the Jacobian of the Lorentz transformation
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables of Weight_Jacobian_Lorentz
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Weight_Jacobian_Lorentz::m_M
    \brief mass of the decaying particle (invariant mass of initial state)
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_Q0
    \brief energy of the charged final states in the \f$ Q_C \f$ rest frame before treatment (\f$ K=0 \f$ event)
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_QN0
    \brief energy of the neutral final states in the \f$ Q_C \f$ rest frame before treatment (\f$ K=0 \f$ event)
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_P0
    \brief energy of the charged final states in the \f$ P_C \f$ rest frame after treatment
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_PN0
    \brief energy of the neutral final states in the \f$ P_C \f$ rest frame after treatment
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_K0
    \brief energy of the generated photons in the \f$ P_C \f$ rest frame
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_mCMSQ
    \brief invariant mass of \f$ Q_C \f$
  */

  /*!
    \var double Weight_Jacobian_Lorentz::m_mCMSP
    \brief invariant mass of \f$ P_C \f$
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods of Weight_Jacobian_Lorentz
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Jacobian_Lorentz::Weight_Jacobian_Lorentz(Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, double, Dipole_Type::code)
    \brief initialises all variables and calculates the weight

    The arguments to be passed are the multipole after treatment, the final state neutral particles after treatment, the multipole before treatment, the final state neutral particles before treatment, the generated photons the invariant mass of the initial state (mass of the decaying particle) and the Dipole_Type of the configuration.

    Calls <tt>Weight_Jacobian_Lorentz::CalculateWeight()</tt> and <tt>Weight_Jacobian_Lorentz::CalculateMax()</tt> to calculate the weight and its maximum. In case of a final state only multipole the latter value is the contribution to the maximum total weight as well as the maximum of \f$ W_{J,L} \f$ by itself, but this is not generally true if the initial state is charged as well.
  */

  /*!
    \fn void Weight_Jacobian_Lorentz::CalculateWeight()
    \brief calculates the weight of the event given

    The weight is \f$ W_{J,L} = \frac{m_C^3}{M^2(P^0+P_N^0+K^0)} \f$.
  */

  /*!
    \fn void Weight_Jacobian_Lorentz::CalculateMax()
    \brief calculates the maximum weight of the configuration

    The maximum weight is calculate at \f$ K=0 \f$, thus \f$ W_{J,L;max} = \frac{m_C^3}{M^2(Q^0+Q_N^0)} \f$.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /*!
    \class Weight_Jacobian_Mapping
    \brief calculates \f$ W_{J,M} \f$, the weight associated with the mapping procedure
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member variables of Weight_Jacobian_Mapping
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Particle_Vector Weight_Jacobian_Mapping::m_newdipole
    \brief contains all charged particles after treatment
  */

  /*!
    \var Particle_Vector Weight_Jacobian_Mapping::m_olddipole
    \brief contains all neutral particles before treatment
  */

  /*!
    \var Particle_Vector Weight_Jacobian_Mapping::m_newspectator
    \brief contains all neutral final state particles after treatment
  */

  /*!
    \var Particle_Vector Weight_Jacobian_Mapping::m_oldspectator
    \brief conatins all neutral final state particles before treatment
  */

  /*!
    \var double Weight_Jacobian_Mapping::m_M
    \brief invariant mass of the initial state (mass of the decaying particle)
  */

  /*!
    \var double Weight_Jacobian_Mapping::m_u
    \brief scaling parameter
  */

  /*!
    \var Vec3D Weight_Jacobian_Mapping::m_K
    \brief total 3-momentum of the photons generated in the \f$ P_C \f$ rest frame
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Weight_Jacobian_Mapping
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Jacobian_Mapping::Weight_Jacobian_Mapping(Particle_Vector, Particle_Vector, Particle_Vector, Particle_Vector, Vec4D, double, double, Dipole_Type::code)
    \brief initialises all variables and calculates the weight

    The arguments to be passed are the multipole after treatment, all final state neutral particles after treatment, the multipole before treatment, all final state neutral particles before treatment, the total momentum of the photons generated, the invariant mass of the initial state, the scaling parameter and the Dipole_Type.

    Calls <tt>Weight_Jacobian_Mapping::CalculateWeight()</tt> and <tt>Weight_Jacobian_Mapping::CalculateMax()</tt> to calculate the weight and its maximum. The latter value is the contribution to the maximum total weight as well as the maximum of \f$ W_{J,M} \f$ by itself.
  */

  /*!
    \fn void Weight_Jacobian_Mapping::CalculateWeight()
    \brief calculates the weight of the event given

    The weight is \f$ W_{J,M} = u^{3n-4}\prod\limits_{i=1}^n\left(\frac{q_i^0}{p_i^0}\right)\frac{\frac{\vec{p}^2}{p^0}-\sum_{C,N}\frac{\vec{q}_i^2}{q_i^0}}{\frac{\vec{p}^\prime\vec{p}}{p^{\prime 0}}-\sum_{C,N}\frac{\vec{p}_i\vec{q}_i}{p_i^0}} \f$.
  */

  /*!
    \fn void Weight_Jacobian_Mapping::CalculateMax()
    \brief calculates the maximum weight of the configuration given

    The maximum weight is 1.
  */
}

#endif
