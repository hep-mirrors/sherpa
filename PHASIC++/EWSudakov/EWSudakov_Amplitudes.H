#ifndef PHASIC_EWSudakov_EWSudakovAmplitudes_H
#define PHASIC_EWSudakov_EWSudakovAmplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"

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

  class Process_Base;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EW Sudakov corrections
   */
  class EWSudakov_Amplitudes {

  public:

    EWSudakov_Amplitudes(Process_Base&);

    ATOOLS::Cluster_Amplitude& Unrotated() noexcept { return *ampls.find(-1)->second; }
    ATOOLS::Cluster_Amplitude& Rotated(size_t legindex);
    std::vector<size_t>& LegPermutation(size_t legindex);

    size_t NumberOfLegs() noexcept { return ampls.find(-1)->second->Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);

  private:

    // helper methods for the construction
    static std::map<int, ATOOLS::Cluster_Amplitude_UP> CreateAmplitudes(Process_Base&);
    static ATOOLS::Cluster_Amplitude_UP CreateAmplitude(Process_Base&);
    static ATOOLS::Cluster_Amplitude_UP CreateSU2RotatedAmplitude(
        const ATOOLS::Cluster_Amplitude_UP&, size_t legindex);
    static std::map<size_t, std::vector<size_t>> CalculateLegPermutations(
        const std::map<int, ATOOLS::Cluster_Amplitude_UP>&);
    static std::vector<size_t> CalculateLegPermutation(
        const ATOOLS::Cluster_Amplitude_UP& ampl);

    /// map flavour-roated leg index -> ampl, -1 gives the unrotated ampl
    std::map<int, ATOOLS::Cluster_Amplitude_UP> ampls;

    /// give permutation of flavours between an ampl and its rotated version
    std::vector<size_t> NewOrdering(const ATOOLS::Cluster_Amplitude_UP&) const;

    /// map flavour-rotated leg index -> permutation
    std::map<size_t, std::vector<size_t>> legpermutations;

  };

}

#endif

