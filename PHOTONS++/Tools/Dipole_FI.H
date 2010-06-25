#ifndef PHOTONS_Tools_Dipole_FI_H
#define PHOTONS_Tools_Dipole_FI_H

#include "PHOTONS++/Tools/Dress_Blob_Base.H"

namespace PHOTONS {
  class Dipole_FI: public Dress_Blob_Base {
    protected:
      double DetermineMaximumPhotonEnergy();

      virtual void   CalculateAvaragePhotonNumber(const double&, const double&);
      virtual bool   CheckIfExceedingPhotonEnergyLimits();
      virtual void   CheckMomentumConservationInQCMS(const ATOOLS::Poincare&,
                                                     const ATOOLS::Poincare&);
      virtual void   CorrectMomenta();
      virtual void   DefineDipole();
      virtual double Func(const double&, const std::vector<double>&,
                          const std::vector<double>&,
                          const std::vector<ATOOLS::Vec3D>&, const double&);
      virtual void   ResetVariables();
      virtual void   ReturnMomenta();

    public:
      Dipole_FI(const Particle_Vector_Vector&);
      virtual ~Dipole_FI();

      virtual void AddRadiation();
  };


  

  /*!
    \file Dipole_FI.H
    \brief contains the class Dipole_FI, the main treatment class for initial-final multipoles
  */

  /*!
    \class Dipole_FI
    \brief main treatment class for initial-final multipoles, daughter of Dress_Blob_Base
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods of Dipole_FI
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Dipole_FI::Dipole_FI(Particle_Vector_Vector)
    \brief initialises main variables necessary for the treatment of initial-final multipoles

    At the moment the implementation is purely for single particle decays. Although 
    most of the treatment is general some points explicitely assume a single charged 
    initial state.
  */

  /*!
    \fn Dipole_FI::~Dipole_FI
    \brief deletes all copies of particles made
  */

  /*!
    \fn double Dipole_FI::DetermineMaximumPhotonEnergy()
    \brief determines the maximum energy of a single photon in the given configuration
  */

  /*!
    \fn void Dipole_FI::AddRadiation()
    \brief manages treatment method, at the end blob is dressed

    The course of action:
    -# call to <tt>Dress_Blob_Base::DefineDipole()</tt> to create the working 
       copies of all particles that may be modified
    -# boost into the \f$ Q_C \f$-system (multipole restframe) with the charged 
       initial state particle in \f$ -z \f$-direction (handy allways, necessary 
       in specialised dipole treatment), definition of \f$ Q,Q_N \f$
    -# calculation of the avarage photon number of the configuration; calls
       <tt>Avarage_Photon_Number</tt> class in case of a multipole or the specialised
       <tt>Dipole_FI::CalculateAvaragePhotonNumber(double,double)</tt> in case 
        of a dipole \n
       \f$ \bar{n}>0 \f$ if \f$ \omega_{max} > \omega_{min} \f$
    -# photon generation procedure (if \f$ \bar{n} > 0 \f$):
      - reset all variables that might have been edited in previous rejected 
        generation cycles to their values before generation
      - generate all photons according to the respective distributions of each property
        (<tt>Dress_Blob_Base::GeneratePhotons</tt>), define \f$ K \f$ and check 
        whether momentum restrictions are fulfilled, i.e. momentum reconstruction 
        is possible (<tt>Dipole_FI::CheckIfExceedingPhotonEnergyLimits()</tt>)
      - if \f$ n > 0 \f$, reconstruct the momenta of all other particles in the 
        new phase space (<tt>Dipole_FI::CorrectMomenta()</tt>), all momenta are now 
        in the \f$ P_C \f$ rest frame
      - calculate the weight of the generated event 
        (<tt>Dress_Blob_Base::CalculateWeights()</tt>) and accept or reject
      - if accepted boost all momenta back into the \f$ Q_C \f$ rest frame
      .
    -# call to <tt>Dipole_FI::CheckMomentumConservationInQCMS()</tt> to check 
       momentum conservation
    -# if \f$ n>0 \f$, boost all momenta back to original frame and paste the 
       modified momenta back into the original particles via 
       <tt>Dipole_FI::ReturnMomenta()</tt>
  */

  /*!
    \fn void Dipole_FI::CalculateAvaragePhotonNumber(double, double)
    \brief calculates the avarage photon number of the given dipole configuration

    This method is only called in case of a dipole. The values to be given are 
    both \f$ \beta \f$ 's in the dipole's rest frame.
  */

  /*!
    \fn bool Dipole_FI::CheckIfExceedingPhotonEnergyLimits()
    \brief checks whether momentum reconstruction is possible
  */

  /*!
    \fn void Dipole_FI::CheckMomentumConservationInQCMS()
    \brief checks momentum conservation in the \f$ Q_C \f$ rest frame, value written to m_success
  */

  /*!
    \fn void Dipole_FI::CorrectMomenta()
    \brief reconstructs all particles' momenta in the new phase space

    The scaling parameter is determined via <tt>Dress_Blob_Base::DetermineU()</tt>. 
    Then all momenta are reconstructed in the \f$ P_C \f$ rest frame.
  */

  /*!
    \fn void Dipole_FI::DefineDipole()
    \brief makes two copies of all charged and final state neutral particles
  */

  /*!
    \fn double Dipole_FI::Func(double, std::vector<double>, std::vector<double>, std::vector<Vec3D>, double)
    \brief returns the value of \f$ f(u) \f$, called by <tt>Dress_Blob_Base::DetermineU()</tt>

    Takes as arguments the squared mass of the initial state particle, the 
    charged final state particles, the neutral final state particles, the 
    3-vectors of all final state particles in the \f$ Q_C \f$ rest frame 
    and the value \f$ u \f$ where \f$ f(u) \f$ shall be evaluated.
  */

  /*!
    \fn void Dipole_FI::ResetVariables()
    \brief resets variables to their values before the correction procedure
  */

  /*!
    \fn void Dipole_FI::ReturnMomenta()
    \brief pastes the new momenta back into the original particles
  */

}

#endif
