#ifndef PHOTONS_Tools_Dress_Blob_Base_H
#define PHOTONS_Tools_Dress_Blob_Base_H

#include "ATOOLS/Math/Vector.H"
#include "PHOTONS++/Main/Dipole_Type.H"

namespace ATOOLS {
  class Poincare;
}

namespace PHOTONS {
  class Dress_Blob_Base {
    protected:
      bool                    m_photonsadded;
      bool                    m_success;
      bool                    m_soft;
      Dipole_Type::code       m_dtype;
      Particle_Vector_Vector  m_pvv;
      Particle_Vector_Vector  m_pvv_new;
      ATOOLS::Particle_Vector m_chargedinparticles;
      ATOOLS::Particle_Vector m_neutralinparticles;
      ATOOLS::Particle_Vector m_chargedoutparticles;
      ATOOLS::Particle_Vector m_neutraloutparticles;
      ATOOLS::Particle_Vector m_softphotons;
      ATOOLS::Particle_Vector m_olddipole;
      ATOOLS::Particle_Vector m_newdipole;
      ATOOLS::Particle_Vector m_oldspectator;
      ATOOLS::Particle_Vector m_newspectator;
      ATOOLS::Particle *      p_newinitialstate;
      ATOOLS::Vec4D           m_K;
      ATOOLS::Vec4D           m_p;
      ATOOLS::Vec4D           m_pN;
      ATOOLS::Vec4D           m_P;
      ATOOLS::Vec4D           m_PN;
      ATOOLS::Vec4D           m_Q;
      ATOOLS::Vec4D           m_QN;
      double                  m_accu;
      double                  m_genweight;
      double                  m_genmaxweight;
      double                  m_nbar;
      int                     m_n;
      double                  m_omegaMax;
      double                  m_omegaMin;
      double                  m_M;
      std::vector<double>     m_mC;
      std::vector<double>     m_mN;
      double                  m_u;

      void                BuildNewParticleVectorVector();
      double              CalculateBeta(const ATOOLS::Vec4D&);
      ATOOLS::Vec4D       CalculateMomentumSum(const ATOOLS::Particle_Vector&);
      void                CalculateWeights();
      void                CheckAvaragePhotonNumberForNumericalErrors();
      void                DeleteAll(ATOOLS::Particle_Vector&);
      void                DeleteAllPhotons();
      void                DetermineU();
      std::vector<double> GenerateNumberAndEnergies();
      void                GeneratePhotons(const double&, const double&);
      void                GeneratePhotons(const IdPairNbarVector&);
      double              KallenFunction(const double&, const double&,
                                         const double&);

      virtual void   CalculateAvaragePhotonNumber(const double&,
                                                  const double&) = 0;
      virtual bool   CheckIfExceedingPhotonEnergyLimits() = 0;
      virtual void   CheckMomentumConservationInQCMS(const ATOOLS::Poincare&,
                                                     const ATOOLS::Poincare&) = 0;
      virtual void   CorrectMomenta() = 0;
      virtual void   DefineDipole() = 0;
      virtual double Func(const double&, const std::vector<double>&,
                          const std::vector<double>&,
                          const std::vector<ATOOLS::Vec3D>&,
                          const double&) = 0;
      virtual void   ResetVariables() = 0;
      virtual void   ReturnMomenta() = 0;

    public:
      Dress_Blob_Base();
      virtual ~Dress_Blob_Base();

      inline bool                     AddedAnything()         { return m_photonsadded; }
      inline bool                     DoneSuccessfully()      { return m_success; }
      inline Particle_Vector_Vector   GetCorrectedParticles() { return m_pvv; }
      inline ATOOLS::Particle *       GetPhoton(int i)        { return m_softphotons.at(i); }
      inline int                      GetPhotonNumber()       { return m_softphotons.size(); }
      inline ATOOLS::Particle_Vector  GetPhotons()            { return m_softphotons; }

      virtual void AddRadiation() = 0;
  };


  

  /*!
    \file Dress_Blob_Base.H
    \brief contains the class Dress_Blob_Base
  */

  /*!
    \class Dress_Blob_Base
    \brief base class for treatment classes Dipole_FF and Dipole_FI
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Dress_Blob_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var bool Dress_Blob_Base::m_photonsadded
    \brief true if more than zero photons were added to the blob, i.e if the momentum configuration has changed
  */

  /*!
    \var bool Dress_Blob_Base::m_success
    \brief true if event generation ended successfully, e.g. momentum conservation is fulfilled
  */

  /*!
    \var bool Dress_Blob_Base::m_soft
    \brief true if only soft corrections should be accounted for

    Set by <tt>Photons::s_mode</tt> which is set by the YFSMODE switch 
    in Fragmentation.dat.
  */

  /*!
    \var Dipole::code Dress_Blob_Base::m_dtype
    \brief specifies the Dipole_Type
  */

  /*!
    \var Particle_Vector_Vector Dress_Blob_Base::m_pvv
    \brief ensemble of particle Vectors passed from Define_Dipole
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_chargedinparticles
    \brief first element of m_pvv, contains all charged initial state particles
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_neutralinparticles
    \brief second element of m_pvv, contains all neutral initial state particles
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_chargedoutparticles
    \brief third element of m_pvv, contains all charged final state particles
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_neutraloutparticles
    \brief fourth element of m_pvv, contains all neutral final state particles
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_softphotons
    \brief contains all generated QED bremsstrahlung

    This Particle_Vector contains all photons generated in the event. If the 
    event is rejected its contents will be deleted and the vector erased. with 
    the routine Dress_Blob_Base::DeleteAll(Particle_Vector). If the event is 
    successful the routine Define_Dipole::AddRadiation() will add its contents 
    to the blob as final state particles.
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_olddipole
    \brief contains a copy of all charged particles of the original blob

    A copy of all charged particles needs to be made since their properties will 
    be modified within the algorithm and the original blob must not be modified 
    in case the event fails.
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_newdipole
    \brief contains a copy of all charged particles after the treatment

    A second copy of all charged particles needs to be made within the generation 
    process since the original configuration must not be modified in case the event 
    is rejected and needs to be regenerated.
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_oldspectator
    \brief contains a copy of all neutral final state particles of the original blob

    A copy of all neutral final state particles needs to be made since their 
    properties will be modified within the algorithm and the original blob must 
    not be modified in case the event fails.
  */

  /*!
    \var Particle_Vector Dress_Blob_Base::m_newspectator
    \brief contains a copy of all neutral final state particles after the treatment

    A second copy of all neutral final state particles needs to be made within 
    the generation process since the original configuration must not be modified 
    in case the event is rejected and needs to be regenerated.
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_K
    \brief sum of all QED bremsstrahlung photon momenta
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_p
    \brief sum of all charged initial state momenta
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_pN
    \brief sum of all neutral initial state momenta
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_P
    \brief sum of all charged final state momenta after radiation
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_PN
    \brief sum of all neutral final state momenta after radiation
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_Q
    \brief sum of all charged final state momenta before radiation
  */

  /*!
    \var Vec4D Dress_Blob_Base::m_QN
    \brief sum of all neutral final state momenta before radiation
  */

  /*!
    \var double Dress_Blob_Base::m_accu
    \brief accuracy used to determine the state of the momentum conservation

    Set by the global accuracy set in ...
  */

  /*!
    \var double Dress_Blob_Base::m_genweight
    \brief weight of the event generated
  */

  /*!
    \var double Dress_Blob_Base::m_genmaxweight
    \brief maximum weight of the configuration for which the event should be generated
  */

  /*!
    \var double Dress_Blob_Base::m_omegaMax
    \brief maximal energy for QED bremsstrahlung, determined kinematically
  */

  /*!
    \var double Dress_Blob_Base::m_omegaMin
    \brief minimal energy for QED bremsstrahlung, infrared cut-off
  */

  /*!
    \var double Dress_Blob_Base::m_nbar
    \brief avarage photon number \f$ \bar{n} \f$
  */

  /*!
    \var double Dress_Blob_Base::m_n
    \brief actual photon number \f$ n \f$ of the event
  */

  /*!
    \var double Dress_Blob_Base::m_M
    \brief mass of the initial state particle, invariant mass of the initial/final state
  */

  /*!
    \var std::vector<double> Dress_Blob_Base::m_mC
    \brief mass of all charged final state particles individually
  */

  /*!
    \var std::vector<double> Dress_Blob_Base::m_mN
    \brief mass of all neutral final state particle individually
  */

  /*!
    \var double Dress_Blob_Base::m_u
    \brief scaling parameter in the momentum reconstruction
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of all member method for Dress_Blob_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Dress_Blob_Base::Dress_Blob_Base
    \brief initialises general event generation variables, e.g infrared cut-off
  */

  /*!
    \fn double Dress_Blob_Base::CalculateBeta(Vec4D)
    \brief calculates \f$ \beta \f$ of a given 4-vector
  */

  /*!
    \fn Vec4D Dress_Blob_Base::CalculateMomentumSum(Particle_Vector)
    \brief calculates the sum of all momenta of the Particle_Vector given
  */

  /*!
    \fn void Dress_Blob_Base::CalculateWeights()
    \brief determines the all weights of the event generated by calling all weight classes

    If <tt>m_soft</tt> is set to true, higher order corrections will 
    not be taken into account.
  */

  /*!
    \fn void Dress_Blob_Base::DeleteAll(Particle_Vector)
    \brief deletes all references of the pointer the given Particle_Vector contains and erases the vector itself
  */

  /*!
    \fn void Dress_Blob_Base::DeleteAllPhotons()
    \brief calls Dress_Blob_Base::DeleteAll(m_softphotons)
  */

  /*!
    \fn void Dress_Blob_Base::DetermineU()
    \brief determines the scaling parameter m_u

    This method determines the scaling parameter \f$ u \f$ as the solution of 
    \f$ f(u) = 0 \f$, which is implemented in the respective treatment classes 
    Dipole_FF and Dipole_FI as <tt>Func(double, std::vector<double>, 
    std::vector<double>, std::vector<Vec3D>, double)</tt>.

    \f$ u \f$ is determined via nested intervals to a precission of \f$ 10^{-12} \f$.
  */

  /*!
    \fn std::vector<double> Dress_Blob_Base::GenerateNumberAndEnergies()
    \brief generates the actual photon number of the event and each photon's energy

    The method first generates the number of photons to be generated according to 
    a Poisson distribution with mean \f$ \bar{n} \f$ and instantaneously generates 
    their energies according to a \f$ \frac{1}{E} \f$ energy distribution. Sets 
    m_n and returns a vector containing the energies.
  */

  /*!
    \fn void Dress_Blob_Base::GeneratePhotons(const double, const double)
    \brief generates a photon configuration in a dipole

    The values to be passed are both \f$ \beta \f$ 's. First calls
     <tt>Dress_Blob_Base::GenerateNumberAndEnergies()</tt> to generate the 
    actual photon number of the event and their energies.  Then it calls the 
    class Generate_One_Photon in its dipole constructor to generate the 
    complete photon and puts it into the m_softphotons Particle_Vector.
  */

  /*!
    \fn void Dress_Blob_Base::GeneratePhotons(std::vector<double>)
    \brief generates a photon configuration in a multipole

    The variable passed is a vector containing the avarage photon numbers 
    \f$ \bar{n}_i \f$ of all individual dipoles (for the ordering convention 
    see Generate_Multipole_Photon_Angle::IndexLookup()). The generation of 
    the actual photon number and their energies proceeds as above, but then 
    instead of the dipole constructor of the Generate_One_Photon class its 
    multipole constructor is called and, again, the complete photons are put 
    into the m_softphotons Particle_Vector.

    The multipole version of GeneratePhotons is able to handle dipoles as well, 
    but the dipole version is faster due to the existence of analytical 
    solutions in this case.
  */

  /*!
    \fn double Dress_Blob_Base::KallenFunction(double, double, double)
    \brief calculates the Kallen-Function of three variables \f$ x,y,z \f$
  */

  /*!
    \fn void Dress_Blob_Base::CalculateAvaragePhotonNumber(const double, const double)
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn bool Dress_Blob_Base::CheckIfExceedingPhotonEnergyLimits()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn void Dress_Blob_Base::CheckMomentumConservationInQCMS()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn void Dress_Blob_Base::CorrectMomenta()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn void Dress_Blob_Base::DefineDipole()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn double Dress_Blob_Base::Func(double, std::vector<double>, std::vector<double>, std::vector<Vec3D>, double)
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn void Dress_Blob_Base::ResetVariables()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn void Dress_Blob_Base::ReturnMomenta()
    \brief purely virtual in Dress_Blob_Base
  */

  /*!
    \fn bool Dress_Blob_Base::AddedAnything()
    \brief returns m_photonsadded
  */

  /*!
    \fn bool Dress_Blob_Base::DoneSuccessfully()
    \brief returns m_success
  */

  /*!
    \fn Particle_Vector_Vector Dress_Blob_Base::GetCorrectedParticles()
    \brief returns m_pvv with the corrected particles (except photons) as entries
  */

  /*!
    \fn Particle * Dress_Blob_Base::GetPhoton(int)
    \brief returns the photon at the given psoition in m_softphotons
  */

  /*!
    \fn int Dress_Blob_Base::GetPhotonNumber()
    \brief returns the number of photon generated in the event
  */

  /*!
    \fn Particle_Vector Dress_Blob_Base::GetPhotons()
    \brief returns m_softphotons
  */

  /*!
    \fn void Dress_Blob_Base::AddRadiation()
    \brief purely virtual in Dress_Blob_Base
  */



}

#endif
