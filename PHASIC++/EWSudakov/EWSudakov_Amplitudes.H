#ifndef PHASIC_EWSudakov_EWSudakovAmplitudes_H
#define PHASIC_EWSudakov_EWSudakovAmplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include <memory>
#include <unordered_set>
#include <unordered_map>

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

  enum class EWSudakov_Amplitude_Type {
    Base,
    ZPhotonInterference,
  };

  class Process_Base;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EW Sudakov corrections
   */
  class EWSudakov_Amplitudes {

    using Cluster_Amplitude_UPM
      = std::unordered_map<size_t, ATOOLS::Cluster_Amplitude_UP>;

  public:

    EWSudakov_Amplitudes(Process_Base*);

    ATOOLS::Cluster_Amplitude& Unrotated() noexcept { return *baseampl; };
    ATOOLS::Cluster_Amplitude& Rotated(EWSudakov_Amplitude_Type, size_t legindex);
    std::vector<size_t>& LegPermutation(EWSudakov_Amplitude_Type, size_t legindex);

    ATOOLS::ClusterAmplitude_Vector AllAmplitudes() noexcept;

    size_t NumberOfLegs() noexcept { return Unrotated().Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);

  private:

    // helper methods for the construction
    Cluster_Amplitude_UPM CreateZPhotonInterferenceAmplitudes(Process_Base*) const;
    static ATOOLS::Cluster_Amplitude_UP CreateAmplitude(Process_Base*);
    static ATOOLS::Cluster_Amplitude_UP CreateZPhotonInterferenceAmplitude(
        const ATOOLS::Cluster_Amplitude_UP&, size_t legindex);
    static std::map<size_t, std::vector<size_t>>
      CalculateLegPermutations(const Cluster_Amplitude_UPM&);
    static std::vector<size_t> CalculateLegPermutation(
        const ATOOLS::Cluster_Amplitude_UP& ampl);

    /// unrotated amplitude
    ATOOLS::Cluster_Amplitude_UP baseampl;

    /// amplitudes for the Z/photon interference
    Cluster_Amplitude_UPM zphotoninterferenceampls;

    /// give permutation of flavours between an ampl and its rotated version
    std::vector<size_t> NewOrdering(const ATOOLS::Cluster_Amplitude_UP&) const;

    /// permutations of legs for each element in zphotoninterferenceampls
    std::map<size_t, std::vector<size_t>> zphotonlegpermutations;

  };

}

#endif

