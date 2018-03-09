#ifndef PHASIC_EWSudakov_EWSudakovAmplitudes_H
#define PHASIC_EWSudakov_EWSudakovAmplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include <array>
#include <memory>
#include <unordered_set>
#include <map>

// TODO: Migrate this to Cluster_Amplitude
namespace ATOOLS {
  auto delAmpl = [](Cluster_Amplitude* ampl) { ampl->Delete(); };
  using Cluster_Amplitude_UP =
    std::unique_ptr<Cluster_Amplitude, decltype(delAmpl)>;
  Cluster_Amplitude_UP MakeClusterAmpl()
  {
    return std::unique_ptr<Cluster_Amplitude, decltype(delAmpl)>{
      Cluster_Amplitude::New(), delAmpl };
  }
  Cluster_Amplitude_UP CopyClusterAmpl(const Cluster_Amplitude_UP& ampl)
  {
    return std::unique_ptr<Cluster_Amplitude, decltype(delAmpl)>{
      ampl->Copy(), delAmpl };
  }
}

namespace PHASIC {

  using Two_Leg_Indizes = std::array<size_t, 2>;

  enum class EWSudakov_Amplitude_Type {
    Base,
    LSCZ,  ///< at leg i, Z is replaced with photon or vice versa
    SSCW   ///< a W moves from i to j, such that i and j are SU(2)-rotated
  };

  class Process_Base;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EW Sudakov corrections
   */
  class EWSudakov_Amplitudes {

    using Leg_Index_Set
      = std::set<size_t>;
    using Cluster_Ampl_Key
      = std::pair<EWSudakov_Amplitude_Type, Leg_Index_Set>;
    using Cluster_Amplitude_UPM
      = std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude_UP>;
    using Permutation_Map
      = std::map<Cluster_Ampl_Key, std::vector<size_t>>;

    static const Cluster_Ampl_Key s_baseamplkey;

  public:

    EWSudakov_Amplitudes(Process_Base*);

    ATOOLS::Cluster_Amplitude& Unrotated() noexcept;
    ATOOLS::Cluster_Amplitude& Rotated(EWSudakov_Amplitude_Type, Leg_Index_Set);
    std::vector<size_t>& LegPermutation(EWSudakov_Amplitude_Type, Leg_Index_Set);

    size_t NumberOfLegs() noexcept { return Unrotated().Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);

    Cluster_Amplitude_UPM::iterator begin() { return ampls.begin(); }
    Cluster_Amplitude_UPM::iterator end() { return ampls.end(); }
    Cluster_Amplitude_UPM::const_iterator cbegin() { return ampls.cbegin(); }
    Cluster_Amplitude_UPM::const_iterator cend() { return ampls.cend(); }

  private:

    // helper methods for the construction
    Cluster_Amplitude_UPM CreateAmplitudes(Process_Base*) const;
    static ATOOLS::Cluster_Amplitude_UP CreateAmplitude(Process_Base*);
    static ATOOLS::Cluster_Amplitude_UP CreateLSCZAmplitude(
        const ATOOLS::Cluster_Amplitude_UP&, size_t legindex);
    static ATOOLS::Cluster_Amplitude_UP CreateSSCWAmplitude(
        const ATOOLS::Cluster_Amplitude_UP&, Two_Leg_Indizes indizes);
    static Permutation_Map CreatePermutations(const Cluster_Amplitude_UPM&);
    static std::vector<size_t> CalculateLegPermutation(
        const ATOOLS::Cluster_Amplitude_UP& ampl);

    /// amplitudes for the Z/photon interference
    Cluster_Amplitude_UPM ampls;

    /// give permutation of flavours between an ampl and its rotated version
    std::vector<size_t> NewOrdering(const ATOOLS::Cluster_Amplitude_UP&) const;

    /// permutations of legs for each element in ampls
    Permutation_Map permutations;

  };

}

#endif
