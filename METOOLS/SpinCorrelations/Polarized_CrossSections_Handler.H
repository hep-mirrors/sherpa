#ifndef METOOLS_SpinCorrelations_Polarized_CrossSections_Handler_H
#define METOOLS_SpinCorrelations_Polarized_CrossSections_Handler_H


#include <vector>
#include <string>
#include <map>
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Spinor.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

namespace ATOOLS{
  class Blob;
  class Particle;
}

namespace METOOLS {
  class Decay_Matrix;
  class Amplitude2_Tensor;
  class PolWeights_Map;
  class Polarization_Vector;
  typedef ATOOLS::Spinor<double> SpinorType;

  class Polarized_CrossSections_Handler {
    std::vector<std::string> m_refsystem;
    std::string m_spinbasis, m_singlepol_channel;
    std::map<std::string, std::string> m_customweights;
    int m_trans_mode;
    ATOOLS::Vec4D m_old_refmom;
    ATOOLS::Vec4D m_new_refmom;
    bool m_helicitybasis;
    bool m_pol_checks;

    void InitPolSettings();
    void InitRefMoms();
    ATOOLS::Vec4D Beta(const ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                       std::string refsystem) const;
    METOOLS::PolWeights_Map* Calculation(ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                                         const std::vector<METOOLS::Decay_Matrix>& decay_matrices,
                                         std::map<int, METOOLS::Polarization_Vector>& default_polarization_vectors,
                                         std::map<int, SpinorType>& default_spinors, ATOOLS::Vec4D beta,
                                         std::string refsystem) const;
    void Tests(const METOOLS::Amplitude2_Tensor* amps, const METOOLS::Amplitude2_Tensor* trafo_pol_amps) const;

  public:
    Polarized_CrossSections_Handler();

    std::vector<METOOLS::PolWeights_Map*> Treat(ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                                                const std::vector<METOOLS::Decay_Matrix>& decay_matrices) const;
    inline std::vector<std::string> GetRefSystems() { return m_refsystem; }
  };


  /*!
      \file Polarized_CrossSections_Handler.H
      \brief Declares the class METOOLS::Polarized_CrossSections_Handler; currently only massive, unstable spin 1
      particles are supported
      */
  /*!
    \class Polarized_CrossSections_Handler
    \brief Class to handle the calculation of polarized cross sections out of an Amplitude2_Tensor and decay matrices,
           currently only unstable, massive spin 1 particles are supported

      This class enables to calculate polarized cross sections out of an Amplitude2_Tensor and decay matrices of the
      particles of interest in various polarization definitions. Currently the ComixDefault spin basis ("ComixDefault",
      the helicity basis ("Helicity")
      and user defined constant spin axes are supported to define the spin basis for massive particles. The polarization
      objects can be calculated in the laboratory ("Lab"), the hard decaying particles center of mass ("COM") and the
      parton-parton frame("PPFr") as well as the restframe of each particle of interest where the spin axis is defined
      in the particle's direction of flight in laboratory frame ("RestFrames"), and the center-of-mass frames defined
      by any combination of hard-process particles. The matrix elements are always taken in the laboratory frame, so all
      polarization objects determined in different frames than that are transformed to the laboratory frame after
      calculation.
      To enable various polarization definitions in one simulation run independent of the matrix
      element generator, the Amplitude2_Tensor and the decay matrices are transformed according to the desired
      polarization definition. It is possible to calculate cross sections for several reference frames in one run.
      The result of the calculation performed here is a PolWeights_Map for each polarization definition containing
      polarization fractions of all polarization definitions. Cross sections are calculated out of that during the event
      output together with all other events weights.
    */
  /*!
  \var std::vector<std::string> m_refsystem
  \brief Vector containing the strings specifing the several reference frames which should be used for polarization
         definition, default: std::vector<std::string>(1, "Lab")
  */
  /*!
  \var std::string m_spinbasis
  \brief String characterizes which spinbasis should be used to calculated the polarized matrix elements, currently
         only one spin basis is supported during one simulation run, default: "Helicity"
  */
  /*!
  \var std::string m_singlepol_channel
  \brief String describes the decay channel which characterizes the only hard decaying particle which should be
         polarized in the desired polarized cross sections, form of the string should be the same as under Channels
         in Hard_Decays scoped setting e.g. 24,-11,12 for W+ decaying into a positron and a neutrino,
         default: "no channel"
  */
  /*!
  \var std::map<std::string, std::string> m_custom_weights
  \brief Map containing the settings for custom weights specified by the user, default: empty map
  */
  /*!
  \var int m_trans_mode
  \brief trans_mode specifies how transverse polarized cross sections should be derived from the "base" polarization
  combinations, 0 = incoherent sum of left and right polarization, 1 = coherent sum of left and right polarization
  including left-right-interference, 2 = cross sections/fractions for both transverse polarization definitions are
  calculated, default: 1
  */
  /*!
  \var ATOOLS::Vec4D m_old_refmom
  \brief Reference momentum used in matrix element generation, defines spin axes for massive particles
  */
  /*!
  \var ATOOLS::Vec4D m_new_refmom
  \brief Reference momentum defining the desired spin axis for massive particles, (0.0, 0.0, 0.0, 0.0) if helicity basis
         is chosen
  */
  /*!
  \var bool m_helicitybasis
  \brief Boolean specifies whether the helicity basis is chosen to define the spin axis
  */
  /*!
  \var bool m_pol_checks
  \brief Boolean specifies whether polarization consistency checks should be done (especially unpol=polsum+int and
   checks whether transformation works properly and is unitary)
  */
  /*!
  \fn Polarized_CrossSections_Handler::Polarized_CrossSections_Handler()
  \brief Constructor of the class

   Constructor of the class, it reads in the polarization settings from the user input and sets the default and the
   desired reference vectors which define the spin axis for massive particles if they are constant; the default
   reference vector is the one used for matrix element calculation with COMIX

  Input:
    - no input
      .
  */
  /*!
  \fn Polarized_CrossSections_Handler::InitPolSettings()
  \brief Helper function for constructor to read in all settings important for the polarization calculation from the user
         input

  Input:
    - no input
      .
  Output:
    - no output
    .
  */
  /*!
  \fn Polarized_CrossSections_Handler::InitRefMoms()
  \brief Helper function for constructor to set the default reference momentum m_old_refmom according to the ComixDefault
         used for the matrix element calculation and the desired one m_new_refmom if it is constant (not the helicity
         basis), the reference momenta fix the spin axis for massive particles
         If further spin bases should be supported, they can be added here by introducing a new keyword and the
         corresponding reference momentum.
  Input:
    - no input
      .
  Output:
    - no output
    .
   */
  /*!
  \fn Polarized_CrossSections_Handler::Treat(ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                                             std::vector<METOOLS::Decay_Matrix> decay_matrices) const
  \brief Method steers the calculation of the polarization fractions with the desired polarization definitions
         This method needs to be called for the actual calculation, it then will call Beta() and Calculate()
  Input:
    - <var>signalblob</var>: (pointer to) the signalblob describing the production of the particles whose polarization
                             should be investigated
    - <var>prod_amps</var>: (pointer to) the Amplitude2_Tensor containing the production matrix elements of the particles
                            of interest, it should contain the amplitudes of all polarization combinations of the matrix
                            element and its complex conjugate.
                            Polarization indexes of particles different than the particles of interest should be summed
                            up.
    - <var>decay_matrices</var>: Vector containing the decay matrices of all particles of interest in case they are
                                 unstable, they are calculated during the spin correlation algorithm in the
                                 Decay_Handler_Base
     .
  Output:
   - vector of pointers to PolWeight_Maps, for each reference system one
   .
  */
  /*!
  \fn Polarized_CrossSections_Handler::Beta(const ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                                            std::string refsystem) const
  \brief Function determines the beta for the Lorentz transformation from the laboratory frame to the desired reference
         frame refsystem

  Beta() returns a zero four-vector if the "RestFrames" frame, the laboratory frame or an undefined frame is
  chosen, the calculated beta depends on momenta, not on particle velocities according to the beta input required in the
  Poincare class;
  if new reference frames with the same beta for all particles should be supported by Sherpa, they can be
  implemented here by introducing a new keyword and the corresponding beta.
  Input:
    - <var>signalblob</var>: (pointer to) the signalblob describing the production of the particles whose polarization
                             should be investigated
    - <var>prod_amps</var>: (pointer to) the Amplitude2_Tensor containing the production matrix elements of the particles
                            of interest, it should contain the amplitudes of all polarization combinations of the matrix
                            element and its complex conjugate.
                            Polarization indexes of particles different than the particles of interest should be summed
                            up.
   - <var>refsystem</var>: specifies the desired reference frame
     .
  Output:
   - beta for the Lorentz transformation from the laboratory frame to the desired reference frame refsystem depending on
     particle momenta
   .
  */
  /*!
  \fn Polarized_CrossSections_Handler::Calculation(ATOOLS::Blob* signalblob, const METOOLS::Amplitude2_Tensor* prod_amps,
                                                   std::vector<METOOLS::Decay_Matrix> decay_matrices,
                                                   std::map<int, METOOLS::Polarization_Vector>& default_polarization_vectors,
                                                   std::map<int, SpinorType>& default_spinors, ATOOLS::Vec4D beta,
                                                   std::string refsystem) const
  \brief Method steers the transformation of the matrix element objects, calculates the actual polarized matrix elements
         is case of intermediate particles and generates a new PolWeights_Map from the results containing the labeled
         polarization fractions

  To support particles different than massive intermediate spin 1 particles the according calculation of the polarization
  objects and the corresponding transformations need to be implemented and then called from here
  Input:
    - <var>signalblob</var>: (pointer to) the signalblob describing the production of the particles whose polarization
                             should be investigated
    - <var>prod_amps</var>: (pointer to) the Amplitude2_Tensor containing the production matrix elements of the particles
                            of interest, it should contain the amplitudes of all polarization combinations of the matrix
                            element and its complex conjugate.
                            Polarization indexes of particles different than the particles of interest should be summed
                            up.
    - <var>decay_matrices</var>: Vector containing the decay matrices of all particles of interest in case they are
                                 unstable, they are calculated during the spin correlation algorithm in the
                                 Decay_Handler_Base
   - <var>default_polarization_vectors</var>: Reference of a map containing the polarization vectors calculated in the
                                              default polarization basis (basis, in which the matrix elements
                                              were calculated) can be empty, in that case the necessary polarization
                                              vectors will be calculated; map key = particle number according to
                                              Sherpa's internal particle numbering
   - <var>default_spinors</var>: Reference of a map containing the spinors calculated in the default polarization basis,
                                 (basis, in which the matrix elements were calculated), can be empty, in that
                                 case, the necessary spinors will be calculated; map key = particle number according to
                                 Sherpa's internal particle numbering
   - <var>beta</var>: beta for the Lorentz transformation from the laboratory frame to the desired reference frame
                      refsystem depending on particle momenta instead of velocities according to the beta required by
                      the Poincare class
   - <var>refsystem</var>: specifies the desired reference frame
     .
  Output:
   - (pointer to) PolWeights_Map
   .
  */
  /*!
  \fn Polarized_CrossSections_Handler::Tests(const METOOLS::Amplitude2_Tensor* amps,
                                             const METOOLS::Amplitude2_Tensor* trafo_pol_amps) const
  \brief Method to test the transformed results
  here
  Input:
    - <var>amps</var>: (pointer to) the Amplitude2_Tensor after contraction with the decay matrices in the polarization
                       basis used for matrix element generation
   - <var>trafo_pol_amps</var>: (pointer to the) Amplitude2_Tensor containing the polarized matrix elements of the
                                particles whose polarization should be investigated, it should contain the amplitudes of
                                all polarization combinations of the matrix element and its complex conjugate.
                                Polarization indexes of particles different than the particles of interest should be
                                summed up. For intermediate particles, trafo_pol_amps should contain the product of
                                production and decay matrix elements
     .
  Output:
   - no output
   .
  */
}

#endif