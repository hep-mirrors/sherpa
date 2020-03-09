#ifndef PHASIC_EWSudakov_EWSudakovAmplitudes_H
#define PHASIC_EWSudakov_EWSudakovAmplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "PHASIC++/EWSudakov/EWSudakov.H"
#include "PHASIC++/Main/Color_Integrator.H"

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
  Cluster_Amplitude_UP CopyClusterAmpl(const Cluster_Amplitude* ampl)
  {
    return std::unique_ptr<Cluster_Amplitude, decltype(delAmpl)>{
      ampl->Copy(), delAmpl };
  }
  Cluster_Amplitude_UP CopyClusterAmpl(const Cluster_Amplitude& ampl)
  {
    return std::unique_ptr<Cluster_Amplitude, decltype(delAmpl)>{
      ampl.Copy(), delAmpl };
  }
}

namespace PHASIC {

  class Process_Base;

  using Cluster_Ampl_Key = Leg_Kfcode_Map;
  using Cluster_Amplitude_UPM =
      std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude_UP>;
  using Cluster_Amplitude_PM =
      std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude*>;
  using Permutation_Map = std::map<Cluster_Ampl_Key, std::vector<size_t>>;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EW Sudakov corrections, when using methods that accept a spin combination,
   * it automatically replaces longitudinal vector bosons with the
   * corresponding Goldstone bosons
   */
  class EWSudakov_Amplitudes {

    static const Cluster_Ampl_Key s_baseamplkey;

  public:

    EWSudakov_Amplitudes(Process_Base*,
                         const std::set<EWSudakov_Log_Type>& _activecoeffs);
    ~EWSudakov_Amplitudes();
    ATOOLS::Cluster_Amplitude& BaseAmplitude() noexcept;
    ATOOLS::Cluster_Amplitude&
    BaseAmplitude(std::vector<int> spincombination);
    Leg_Kfcode_Map GoldstoneBosonReplacements(std::vector<int> spincombination);
    ATOOLS::Cluster_Amplitude& SU2TransformedAmplitude(const Leg_Kfcode_Map&);
    std::vector<size_t>& LegPermutation(const Leg_Kfcode_Map&);

    size_t NumberOfLegs() noexcept { return BaseAmplitude().Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);
    void UpdateColors(const Int_Vector& I, const Int_Vector& J);
    double MandelstamS();
    double MandelstamT();
    double MandelstamU();

    const Cluster_Amplitude_PM& All() const { return all_ampls; }
    const Cluster_Amplitude_PM& GoldstoneOnly() const
    {
      return goldstone_only_ampls;
    }

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

    /// all amplitudes that might be needed for EW Sudakov corrections
    Cluster_Amplitude_UPM ampls;
    Cluster_Amplitude_PM all_ampls, goldstone_only_ampls;

    /// give permutation of flavours between an ampl and its transformed version
    std::vector<size_t> NewOrdering(const ATOOLS::Cluster_Amplitude_UP&) const;

    /// permutations of legs for each element in ampls
    Permutation_Map permutations;

    /// used to put all momenta on-shell for transformed matrix elements
    ATOOLS::Momenta_Stretcher stretcher;

    /// this container is only used to accommodate the Momenta_Stretcher API
    ATOOLS::Particle_Vector particles;

  };

}

#endif
