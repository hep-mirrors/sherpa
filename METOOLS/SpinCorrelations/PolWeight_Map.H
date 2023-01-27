#ifndef METOOLS_SpinCorrelations_PolWeight_Map_H
#define METOOLS_SpinCorrelations_PolWeight_Map_H

#include <map>
#include <set>
#include <vector>
#include <memory>

#include "ATOOLS/Math/MyComplex.H"

namespace METOOLS {
  class Amplitude2_Tensor;

  class PolWeight_Map : public std::map<std::string, Complex> {
    Complex m_unpolcrosssec;
    bool m_massive_vb;
    std::map<std::string, std::string> m_custom_weights;
    std::string m_singlepol_channel;
  public:
    PolWeight_Map();
    PolWeight_Map(const METOOLS::Amplitude2_Tensor* amps,
                  std::map<std::string, std::string> custom_weights=std::map<std::string, std::string>(),
                  std::string singlepol_channel="no channel");
    ~PolWeight_Map();
    void Init(const METOOLS::Amplitude2_Tensor* amps, std::string mode="start", std::string spin_label="");
    std::set<std::string> ListofKeys() const;
    PolWeight_Map Transverse (int level) const;
    PolWeight_Map AddCustomWeights (const METOOLS::Amplitude2_Tensor* amps) const;
    PolWeight_Map Unpol(const METOOLS::Amplitude2_Tensor* amps, std::vector<int> particle_numbers) const;
    PolWeight_Map AddSinglePolWeights(const METOOLS::Amplitude2_Tensor* amps) const;
  };

  /*!
    \file PolWeight_Map.H
    \brief Declares the class METOOLS::PolWeight_Map; all methods of this class can be used for spin 1 and spin 1/2
           particles; class inherits from std::map
    */
  /*!
    \class PolWeight_Map
    \brief Class to extract Polarization weights from an Amplitude2_Tensor

      This class enables to extract the polarization weights inclusive interference from an
      Amplitude2_Tensor amps, assgins a reasonable label to every weight which can be used for storing them in a
      WeightMap later and save them in a map
    */
  /*!
\var Complex m_unpolcrosssec
\brief Unpolarized cross section of the Amplitude2_Tensor for which the PolWeightMap is calculated
*/
  /*!
  \var bool m_massive_vb
  \brief Boolean specifies whether VB are included in Amplitude2_Tensor, to decide later, whether transversal
         weights should be calculated or not
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

  \fn PolWeight_Map::PolWeightsMap()
  \brief Constructor of the class, defining an empty map
  Input:
    - no input
      .
  */
  /*!
  \fn PolWeight_Map::PolWeightsMap(METOOLS::Amplitude2_Tensor* amps,
   std::map<std::string, std::string> custom_weights=std::map<std::string, std::string>(),
   std::string singlepol_channel="no channel")
  \brief Constructor of the class

   Constructor of the class generating the PolWeight_Map from given Amplitude2_Tensor by help of Helper function Init
   where also custom weights specified by m_custom_weights are added to the map

  Input:
    - <var>amps</var>: (pointer to an) amplitude2_tensor containing all polarized matrixelements of interest
    - <var>custom_weights</var>: contains weights specified by the user stored in a map with key=Setting name in
                                 YAML-File (Weight, Weightn), value=user specified weight (comma-separated weight
                                 names or particle numbers)
    - <var>singlepol_channel</var>: String describes the decay channel which characterizes the only hard decaying
                                    particle which should be polarized in the desired polarized cross sections,
                                    form of the string should be the same as under Channels
                                    in Hard_Decays scoped setting e.g. 24,-11,12 for W+ decaying into a positron and a
                                    neutrino
      .
  */
  /*!
  \fn PolWeight_Map::Init(METOOLS::Amplitude2_Tensor* amps, std::string mode="start", std::string spin_label="")
  \brief Helper function for constructor to generate basic PolWeightMap

  Extracts all polarized matrixelements from the provided Amplitude2_Tensor amps, normalizes them with the
  unpolarized result (sum over all tensor elements), provides them with a label and saves them in a map (key = label,
  value = polarization fraction); all interference contributions are identified, summed and additionally added to the
  generated map describing the respective polarization combination (+ = positive helicity, - = negative helicity,
  0 = longitudinal); additionally some consistency checks are directly implemented
  Input:
    - <var>amps</var>: (pointer to an) amplitude2_tensor containing all polarized matrixelements of interest
    - optional <var>mode</var>: mode describing whether the function is called the first time (consistency checks are
                                only done in that case) or whether the current fraction is already identified as
                                interference
                                default: "start"
    - optional <var>spin_label</var>: the label for the current polarization fraction is build recursively, during
                                      the single steps of recursion the already generated label is saved in and
                                      propagated by the help of this variable
                                      default: empty string ""
      .
  Output:
    - no output
    .
  */
  /*!
  \fn PolWeight_Map::ListofKeys()
  \brief Generate a list of all keys in a polarization map

  Generate a list of all keys in a polarization map since C++ standard library map does not provide such a
       functionality
  Input:
    - no input
      .
  Output:
    - List of keys as std::vector<std::string>
    .
  */
  /*!
  \fn PolWeight_Map::Transverse(int level)
  \brief Calculate transversal cross sections if vector boson polarizations are considered by summing contributions
         of +- and --helicity together to one additional weight
  Input:
  - <var>level</var>: Method goes through the spin labels by recursive calls of itself, level describes which
                      particle in the label the current method call is working on
  .
  Output:
  - PolWeight_Map containing all transversal weights (the summation of + and - helicity is done for all VB here, so
    there are no weights described by a mixture between transversal and +, - contributions from different VB, for
    that the custom weight setting could be used)
  .
  */
  /*!
  \fn PolWeight_Map::AddCustomWeights(METOOLS::Amplitude2_Tensor* amps)
  \brief Calculate user specified additional weights if weight names are given, otherwise it calls the Unpol method to
         calculate weights where one or several particles are considered as unpolarized
  - <var>amps</var>: (pointer to an) amplitude2_tensor containing all polarized matrixelements of interest
  .
  Output:
  - PolWeight_Map containing all custom weights
  .
  */
  /*!
  \fn PolWeight_Map::Unpol(METOOLS::Amplitude2_Tensor* amps, std::vector<int> particle_numbers,
                           PolWeight_Map temp_map=PolWeight_Map())
  \brief Calculate user specified additional weights if weight names are given, otherwise it calls the Unpol method to
         calculate weights where one or several particles are considered as unpolarized
  - <var>amps</var>: (pointer to an) amplitude2_tensor containing all polarized matrixelements of interest
  - <var>particle_numbers</var>: vector containing the number of that particles which should be considered as
                                 unpolarized; numbers according to particle numbering in YAML
  .
  Output:
  - PolWeight_Map containing all desired polarization weights with partially unpolarized particles
  .
  */
  /*!
  \fn PolWeight_Map::AddSinglePolWeights(METOOLS::Amplitude2_Tensor* amps)
  \brief Calculate single polarized cross sections, single polarized particle is specified by its decay channel
  - <var>amps</var>: (pointer to an) amplitude2_tensor containing all polarized matrixelements of interest
  .
  Output:
  - PolWeight_Map containing all desired single polarization weights
  .
  */
}

#endif
