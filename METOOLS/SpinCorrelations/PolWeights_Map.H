#ifndef METOOLS_SpinCorrelations_PolWeight_Map_H
#define METOOLS_SpinCorrelations_PolWeight_Map_H

#include <map>
#include <set>
#include <vector>
#include <memory>

#include "ATOOLS/Math/MyComplex.H"

namespace METOOLS {
  class Amplitude2_Tensor;

  class PolWeights_Map : public std::map<std::string, Complex> {
    Complex m_unpolcrosssec;
    bool m_massive_vb, m_pol_checks;
    int m_trans_mode;
    std::map<std::string, std::string> m_custom_weights;
    std::string m_singlepol_channel;
    std::set<std::string> m_interference_weights;
    PolWeights_Map* p_all_weights;
  public:
    PolWeights_Map();
    PolWeights_Map(const METOOLS::Amplitude2_Tensor* amps, int trans_mode,
                   std::map<std::string, std::string> custom_weights=std::map<std::string, std::string>(),
                   std::string singlepol_channel="no channel", bool pol_checks=false);
    ~PolWeights_Map();

    std::set<std::string> ListofKeys() const;
    std::string ShortName(std::string name) const;
    std::set<std::string> TransverseKeys(std::set<std::string> keys, int level, bool coherent) const;
    std::vector<std::string> ExpandLabels(const std::vector<std::string>& transverse_labels, int level, int num_particles=0) const;

    void Tests(std::string prefix="");

    void LabelAndSeparate(const METOOLS::Amplitude2_Tensor* amps, const std::string& mode="start",
                          const std::string& prefix="", bool nonzero_weight= true, std::string spin_label="");
    void Transverse (bool coherent);
    void AddCustomWeights (const METOOLS::Amplitude2_Tensor* amps,
                           const std::vector<std::string>& finished_custom_weights);
    std::vector<std::string> Unpol(const METOOLS::Amplitude2_Tensor* amps,
                                   const std::vector<int>& unpol_particle_numbers= std::vector<int>(),
                                   bool non_zero_weight=true);
    void AddSinglePolWeights(const METOOLS::Amplitude2_Tensor* amps);
  };

  /*!
    \file PolWeights_Map.H
    \brief Declares the class METOOLS::PolWeights_Map; all methods of this class can be used for spin 1 and spin 1/2
           particles; class inherits from std::map
    */
  /*!
    \class PolWeights_Map
    \brief Class to extract polarization fractions from an Amplitude2_Tensor

      This class enables to extract polarization fractions and interference contributions from an
      Amplitude2_Tensor amps, assigns a reasonable label to every of those weights and save them in a C++ map;
      beside fractions corresponding to each possible polarization combination and a totaled interference contribution
      (= base weights) it can also calculate transverse polarization fractions if massive VBs are involved as
      intermediate particles and user-specified polarization fractions (partially unpolarized or sum of
      several base weights)
    */
  /*!
  \var Complex m_unpolcrosssec
  \brief Unpolarized amplitude resulting from the sum over all entries of Amplitude2_Tensor for which the PolWeights_Map
         is calculated
  */
  /*!
  \var bool m_massive_vb
  \brief Specifies whether massive VBs are included in Amplitude2_Tensor, if true, transverse weights are calculated
  */
  /*!
  \var int m_trans_mode
  \brief Trans_mode specifies how transverse polarized cross sections should be derived from the "base" polarization
  combinations, 0 = incoherent sum of left and right polarization, 1 = coherent sum of left and right polarization
  including left-right-interference terms, 2 = cross sections/fractions for both transverse polarization definitions are
  calculated
  default: 1
  */
  /*!
  \var bool m_pol_checks
  \brief Boolean specifies whether polarization consistency checks should be done (especially unpol=polsum+int and
   checks whether transformation works properly and is unitary)
  */
  /*!
  \var std::map<std::string, std::string> m_custom_weights
  \brief Map containing the settings for custom weights specified by the user
  */
  /*!
  \var std::string m_singlepol_channel
  \brief String describes the decay channel which characterizes the only hard decaying particle which should be
         polarized in the desired polarized cross sections, form of the string should be the same as under Channels
         in Hard_Decays scoped setting e.g. 24,-11,12 for W+ decaying into a positron and a neutrino
  */
  /*!
  \var std::set<std::string> m_interference_weights
  \brief Set contains all spin labels which describe the interference terms
  */
  /*!
  \var PolWeights_Map* p_all_weights
  \brief Pointer to a PolWeights_Map containing all fractions corresponding to all entries of the Amplitude2_Tensor and
   partially unpolarized fractions derived from that, this includes especially all interference terms
  */
  /*!

  \fn PolWeights_Map::PolWeightsMap()
  \brief Constructor of the class, defining an empty map
  Input:
    - no input
      .
  */
  /*!
  \fn PolWeights_Map::PolWeightsMap(METOOLS::Amplitude2_Tensor* amps,
   std::map<std::string, std::string> custom_weights=std::map<std::string, std::string>(),
   std::string singlepol_channel="no channel", bool pol_checks=false)
  \brief Constructor of the class

   Constructor of the class generating the PolWeights_Map from given Amplitude2_Tensor, from the base weights
   (=normalized entries of amps) also transverse (if massive VBs are involved as intermediate particles) and custom
   weights (if desired) are calculated and added to the PolWeightsMap
   Input:
    - <var>amps</var>: (pointer to an) Amplitude2_tensor containing all polarized matrix elements of interest
    - <var>custom_weights</var>: contains custom weights settings specified by the user, they are stored in a map with
                                 key=Setting name in YAML-File (Weight, Weightn), value=user specified weight
                                 (comma-separated weight names or particle numbers)
    - <var>singlepol_channel</var>: String describes the decay channel which characterizes the only hard decaying
                                    particle which should be polarized in the desired polarized cross sections,
                                    form of the string should be the same as under Channels in Hard_Decays scoped
                                    setting e.g. 24,-11,12 for W+ decaying into a positron and a neutrino
    - <var>pol_checks</var>: Boolean specifies whether polarization consistency checks should be done (especially
                             unpol=polsum+int and checks whether transformation works properly and is unitary)
                             default: false
      .
  */
  /*!
   \fn PolWeights_Map::ListofKeys()
  \brief Generate a list of all keys in a PolWeights_Map
  Input:
    - no input
      .
  Output:
    - std::set<std::string>: List of keys
    .
  */
  /*!
  \fn   PolWeights_Map::ShortName(std::string name)
  \brief Reduces all doubly appearing polarization indices in spin labels to one single index
  Input:
    - <var>name</var>: spin label which should be shortened      .
  Output:
    - std::string: reduced spin label
    .
  */
  /*!
  \fn PolWeights_Map::TransverseKeys(int level, bool coherent, std::set<std::string> keys)
  \brief Generate spin labels describing all possible polarization combinations which include at least one transverse
         polarized, massive VB from base spin labels where the transverse polarizations are separated in left(-)- and
         right-handed(+) polarization
  Input:
  - <var>keys</var>: Set of basic spin labels, the method should use to determine the spin labels where left- and right-
                     handed polarization of massive VBs are summed up; during recursive calls of the method, keys
                     contain the intermediate spin labels which are calculated in the call of the method beforehand
                     (can still contain left-handed- and right-handed-polarization indices for VBs at levels > current
                     level)
  - <var>level</var>: Method goes through the spin labels by recursive calls of itself, level describes which
                      particle in the label the current method call is working on
  - <var>coherent</var>: specifies whether spin labels should describe transverse polarized particles where the
                         corresponding transverse signal includes left-right interference (coherent) or not
                         (incoherent), they are distinguished by using capital (coherent) and small letter (incoherent)
                         T,t respectively                     .
  Output:
  - std::set<std::string>: Set of spin labels describing all possible polarization combinations which includes at least
    one transverse polarized
  .
  */
  /*!
  \fn PolWeights_Map::ExpandLabels(const std::vector<std::string>& strings, int level, int num_particles=0)
  \brief Method returns the spin labels of all polarization fractions which need to be added to receive the polarization
         fraction which corresponds to the passed transverse spin label;
         this is achieved by replacing recursively each T or t in the spin label of interest by contributing + and -
         polarization combinations (++, --  for incoherent definition, ++, --, +-,-+ for coherent definition, first
         index refers to the polarization in the matrix element, the second in its complex conjugate
  Input:
  - <var>strings</var>: Vector of spin labels where at least one massive VB is transversely polarized; if the method
                        is called (initially), this vector should only contain one spin label; the vector object is only
                        chosen here for the recursive calls of the method, during which it is filled with the spin
                        labels to add up
  - <var>level</var>: Method goes through the spin labels by recursive calls of itself, level describes which
                      particle in the label the current method call is working on
  - optional <var>num_particles</var>: Number of (intermediate) particles in the spin label; only for recursive calls,
                                       is determined during the initial call of the method
                                       default: 0
  Output:
  - std::vector<std::string>: Vector of spin labels describing all polarization fractions which need to be added to
    receive the polarization fraction which corresponds to the passed transverse spin label
  .
  */
  /*!
  \fn PolWeights_Map::Tests(std::string prefix=""))
  \brief Method performs some checks after the polarization weights are labeled and added to the PolWeights_Map
  Input:
  - <var>prefix</var>: String specifies possible prefix in the weight name (e.g. "dc" or "Weight1")
                       default: empty string
  Output:
  - no output
    .
  */
  /*!
  \fn PolWeights_Map::LabelAndSeparate(const METOOLS::Amplitude2_Tensor* amps, const std::string& mode="start",
   const std::string& prefix="", bool nonzero_weight= true, std::string spin_label="")
  \brief Method to label the entries of an Amplitude2_Tensor, to calculate polarization fractions from that and
         to separate contributions where all particles are in a definite polarization state (= polarized contributions)
         from terms describing interferences between different polarizations

  Extracts all polarized matrix elements from the provided Amplitude2_Tensor amps, normalizes them with the
  unpolarized result (sum over all Amplitude2_Tensor elements), provides them with a label and add them to the
  PolWeights_Map p_all_weights is pointing to and in case of the polarized contributions also to the PolWeights_Map
  this method is called on (key = label, value = polarization fraction);
  the idea behind this procedure is that the actual constructed PolWeights_Map only contains weights which should be
  printed out in the end while p_all_weights is necessary to include some interference terms into the coherent
  transverse signal definition;
  furthermore, all interference contributions are identified, summed and additionally added to the PolWeights_Map,
  this method is called on;
  general form of the spin labels:
    <particle1>.<polarization_index_in_matrix_element><polarization_index_in_complex_conjugate_matrix_element>_
    <particle2>.<polarization_index_in_matrix_element><polarization_index_in_complex_conjugate_matrix_element> ...
  order of the intermediate particles in the label is according to the order of the particles in the Amplitude2_Tensor
  Input:
    - <var>amps</var>: (pointer to an) Amplitude2_Tensor containing all polarized matrix elements of interest
    - optional <var>mode</var>: mode describing whether the function is called the first time or whether the current
                                fraction is already identified as interference
                                default: "start"
    - optional <var>prefix</var>: specifies a string which should be added at the beginning of each spin label generated
                                  by this method
                                  default: empty string
    - optional <var>nonzero_weight</var>: specifies whether the fractions corresponding to the generated spin labels
                                          should be set to zero (=false), only necessary if partially unpolarized
                                          cross sections for identical particles should be calculated where the
                                          particles considered as polarized are characterized by a certain decay channel
                                          (see AddSinglePolWeights method for details)
                                          default: true
    - optional <var>spin_label</var>: the label for the current polarization fraction is build recursively, during
                                      the single steps of recursion the already generated label is saved in and
                                      propagated by the help of this variable
                                      default: empty string ""
  Output:
    - no output
    .
  */
  /*!
  \fn PolWeights_Map::Transverse(bool coherent)
  \brief Calculate transverse polarization fractions, if VB polarizations are considered, this is done
         by summing contributions of left- and right-handed polarization as well as left-right-interference
         (for coherent definition) together to one additional weight
  Input:
  - <var>coherent</var>: boolean specifies which definition of the transverse polarization should be considered
                         (true: coherent definition including left-right interference, false: incoherent definition, no
                         interference included)
  .
  Output:
  - no output
  .
  */
  /*!
  \fn PolWeights_Map::AddCustomWeights(const METOOLS::Amplitude2_Tensor* amps,
                           const std::vector<std::string>& finished_custom_weights)
  \brief Calculate user specified additional polarization fractions if spin labels are given;
         polarization fractions that correspond to the user specified spin labels are totaled and added as additional
         weight to the PolWeights_Map
  - <var>amps</var>: (pointer to an) Amplitude2_Tensor containing all polarized matrix elements of interest
  - <var>finished_custom_weights</var>: Custom weights which are already calculated beforehand (e.g. partially
                                        unpolarized weights), necessary to test, whether each user input regarding
                                        customized polarized cross sections is valid
  .
  Output:
  - no output
  .
  */
  /*!
  \fn PolWeights_Map::Unpol(const METOOLS::Amplitude2_Tensor* amps,
                            const std::vector<int>& unpol_particle_numbers= std::vector<int>(),
                            bool non_zero_weight=true)
  \brief Calculate polarization fractions where one or several intermediate particles are considered as unpolarized
  - <var>amps</var>: (pointer to an) Amplitude2_Tensor containing all polarized matrix elements of interest
  - <var>unpol_particle_numbers</var>: vector containing the number of those particles which should be considered as
                                       unpolarized; numbers according to particle numbering in Sherpa;
                                       unpol_particle_numbers being not empty is only necessary, if m_custom_weights
                                       does not contain the particle numbers of those particles which should be
                                       considered as unpolarized, e.g. if partially unpolarized weights should be
                                       calculated for processes with identical particles involved and particles
                                       decaying via a certain decay channel should be considered as polarized
  .
  Output:
  - std::vector<std::string> Vector of user inputs specifying which custom weights of the desired ones are already
    calculated by this method
  .
  */
  /*!
  \fn PolWeights_Map::AddSinglePolWeights(METOOLS::Amplitude2_Tensor* amps)
  \brief Calculate single polarization fractions for processes with identical intermediate particles like in
   same sign VB pair production, single polarized particle is specified by its decay channel
   private version, no yet included in manual
  - <var>amps</var>: (pointer to an) Amplitude2_Tensor containing all polarized matrix elements of interest
  .
  Output:
  - PolWeights_Map containing all desired single polarization weights
  .
  */
}

#endif
