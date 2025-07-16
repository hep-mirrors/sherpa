#ifndef EWSudakov_Amplitudes_H
#define EWSudakov_Amplitudes_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "AddOns/EWSud/EWSud.H"
#include "PHASIC++/Main/Color_Integrator.H"

#include <memory>
#include <map>

namespace PHASIC {
  class Process_Base;
}

namespace EWSud {

  using Cluster_Ampl_Key = Leg_Kfcode_Map;
  using Cluster_Amplitude_UPM =
      std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude_UP>;
  using Cluster_Amplitude_PM =
      std::map<Cluster_Ampl_Key, ATOOLS::Cluster_Amplitude*>;

  /**
   * helper class that generates and stores the amplitudes needed to calculate
   * EWSud corrections, when using methods that accept a spin combination,
   * it automatically replaces longitudinal vector bosons with the
   * corresponding Goldstone bosons
   */
  class Amplitudes {

    static const Cluster_Ampl_Key s_baseamplkey;

  public:
    Amplitudes(PHASIC::Process_Base*,
               const std::set<EWSudakov_Log_Type>& _activecoeffs);
    ~Amplitudes();
    ATOOLS::Cluster_Amplitude& BaseAmplitude() noexcept;
    ATOOLS::Cluster_Amplitude&
    BaseAmplitude(std::vector<int> spincombination);
    Leg_Kfcode_Map GoldstoneBosonReplacements(std::vector<int> spincombination);
    ATOOLS::Cluster_Amplitude& SU2TransformedAmplitude(const Leg_Kfcode_Map&);

    size_t NumberOfLegs() noexcept { return BaseAmplitude().Legs().size(); }

    void UpdateMomenta(const ATOOLS::Vec4D_Vector& mom);
    void UpdateColors(const PHASIC::Int_Vector& I, const PHASIC::Int_Vector& J);
    double MandelstamS();
    double MandelstamT();
    double MandelstamU();

    const Cluster_Amplitude_PM& All() const { return all_ampls; }
    const Cluster_Amplitude_PM& GoldstoneOnly() const
    {
      return goldstone_only_ampls;
    }

    // helper methods for the construction
    Cluster_Amplitude_UPM CreateAmplitudes(
        PHASIC::Process_Base*,
        const std::set<EWSudakov_Log_Type>& activecoeffs) const;
    static ATOOLS::Cluster_Amplitude_UP CreateAmplitude(PHASIC::Process_Base*);
    static ATOOLS::Cluster_Amplitude_UP
    CreateSU2TransformedAmplitude(const ATOOLS::Cluster_Amplitude_UP& ampl,
                                  const Leg_Kfcode_Map_Signed&);

    size_t StretcherFailFlag = (1<<4);

  private:

    /// all amplitudes that might be needed for EWSud corrections
    Cluster_Amplitude_UPM ampls;
    Cluster_Amplitude_PM all_ampls, goldstone_only_ampls;

    /// used to put all momenta on-shell for transformed matrix elements
    ATOOLS::Momenta_Stretcher stretcher;

    /// this container is only used to accommodate the Momenta_Stretcher API
    ATOOLS::Particle_Vector particles;

  };

}

#endif
