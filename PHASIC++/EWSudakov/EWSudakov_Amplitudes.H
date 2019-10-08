#ifndef PHASIC_EWSudakov_EWSudakovAmplitudes_H
#define PHASIC_EWSudakov_EWSudakovAmplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/EWSudakov/EWSudakov.H"

#include <memory>
#include <map>

// TODO: Migrate this to Cluster_Amplitude on master and use it everywhere
// when doing this we can get rid of custom smart pointers on master
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

  class Process_Base;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EW Sudakov corrections, when using methods that accept a spin combination,
   * it automatically replaces longitudinal vector bosons with the
   * corresponding Goldstone bosons
   */
  class EWSudakov_Amplitudes {

    using Cluster_Ampl_Key = Leg_Kfcode_Map;
    using Cluster_Amplitude_UPM
      = std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude_UP>;
    using Permutation_Map
      = std::map<Cluster_Ampl_Key, std::vector<size_t>>;

    static const Cluster_Ampl_Key s_baseamplkey;

  public:

    EWSudakov_Amplitudes(Process_Base*,
                         const std::set<EWSudakov_Log_Type>& _activecoeffs);

    ATOOLS::Cluster_Amplitude& BaseAmplitude() noexcept;
    ATOOLS::Cluster_Amplitude&
    BaseAmplitude(std::vector<int> spincombination) noexcept;
    ATOOLS::Cluster_Amplitude& SU2TransformedAmplitude(const Leg_Kfcode_Map&);
    std::vector<size_t>& LegPermutation(const Leg_Kfcode_Map&);

    size_t NumberOfLegs() noexcept { return BaseAmplitude().Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);
    double MandelstamS();
    double MandelstamT();
    double MandelstamU();

    Cluster_Amplitude_UPM::iterator begin() { return ampls.begin(); }
    Cluster_Amplitude_UPM::iterator end() { return ampls.end(); }
    Cluster_Amplitude_UPM::const_iterator cbegin() { return ampls.cbegin(); }
    Cluster_Amplitude_UPM::const_iterator cend() { return ampls.cend(); }

  private:

    // helper methods for the construction
    Cluster_Amplitude_UPM CreateAmplitudes(
        Process_Base*, const std::set<EWSudakov_Log_Type>& activecoeffs) const;
    static ATOOLS::Cluster_Amplitude_UP CreateAmplitude(Process_Base*);
    static ATOOLS::Cluster_Amplitude_UP
    CreateSU2TransformedAmplitude(const ATOOLS::Cluster_Amplitude_UP& ampl,
                                  const Leg_Kfcode_Map_Signed&);
    static Permutation_Map CreatePermutations(
        const Cluster_Amplitude_UPM&);
    static std::vector<size_t> CalculateLegPermutation(
        const ATOOLS::Cluster_Amplitude_UP& ampl);

    /// amplitudes for the Z/photon interference
    Cluster_Amplitude_UPM ampls;

    /// give permutation of flavours between an ampl and its transformed version
    std::vector<size_t> NewOrdering(const ATOOLS::Cluster_Amplitude_UP&) const;

    /// permutations of legs for each element in ampls
    Permutation_Map permutations;

  };

}

#endif
