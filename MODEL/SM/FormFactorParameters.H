#ifndef METOOLS_FormFactorParameters_H
#define METOOLS_FormFactorParameters_H

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"

namespace METOOLS {

  // ------------------------------------------------------------------
  // Resonance type enum – covers all lineshapes we might use
  // ------------------------------------------------------------------
  enum class ResonanceType {
    FixedBreitWigner, // constant width
    GounarisSakurai,  // GS running width (for rho → ππ)
    RunningWidth,     // generic running width (not yet implemented)
    Custom            // user-defined (e.g., for J/ψ, φ, etc.)
  };

  // ------------------------------------------------------------------
  // A single resonance parameter set: mass, width, amplitude, phase,
  // and the type of propagator to use.
  // ------------------------------------------------------------------
  class ResonanceParameters {
  private:

    std::string m_block; // settings block (e.g., "PION_FORM_FACTOR")
    std::string m_sub;   // optional sub-block (e.g., "rho(770)")
    std::string m_tag;   // resonance name (for output)
    double m_m0, m_G0, m_c0, m_phi0; // defaults
    double m_m, m_G, m_c, m_phi;     // current values
    ResonanceType m_type;            // propagator type
    bool m_shown;                    // already printed?
    kf_code m_flav;                  // flavour (for lineshape lookup)

  public:

    // Constructor for simple (no sub-block)
    ResonanceParameters(const std::string& block, const std::string& tag,
                        double m, double Gamma, double c, double phi,
                        ResonanceType type = ResonanceType::FixedBreitWigner,
                        kf_code flav = kf_none);

    // Constructor with sub-block
    ResonanceParameters(const std::string& block, const std::string& sub,
                        const std::string& tag, double m, double Gamma,
                        double c, double phi,
                        ResonanceType type = ResonanceType::FixedBreitWigner,
                        kf_code flav = kf_none);

    // Read from Settings (idempotent)
    void Read();

    // Accessors
    inline double Mass() const { return m_m; }
    inline double Width() const { return m_G; }
    inline Complex Weight() const
    {
      return m_c * Complex(std::cos(m_phi), std::sin(m_phi));
    }
    inline ResonanceType Type() const { return m_type; }
    inline kf_code Flavour() const { return m_flav; }
    inline const std::string& Tag() const { return m_tag; }

    // For printing
    void PrintLine() const;
  };

  // ------------------------------------------------------------------
  // Parameter registry – holds all resonance parameters, grouped by
  // sector (e.g., "Pion", "Kaon", "ThreePion", "TwoPhoton", "Baryon").
  // ------------------------------------------------------------------
  class ParameterRegistry {
  private:

    std::map<std::string,std::vector<ResonanceParameters*> > m_data;

    // Singleton
    ParameterRegistry() = default;
    ~ParameterRegistry() = default;
    ParameterRegistry(const ParameterRegistry&) = delete;
    ParameterRegistry& operator=(const ParameterRegistry&) = delete;

  public:

    static ParameterRegistry& Instance();

    // Add parameters to a sector
    void Add(const std::string& sector, ResonanceParameters* params);

    // Get all parameters for a sector
    const std::vector<ResonanceParameters*>&
    Get(const std::string& sector) const;

    // Read all parameters in a sector (calls Read() on each)
    void ReadSector(const std::string& sector);

    // Print a sector header/table
    void PrintSector(const std::string& sector, const std::string& title = "");
  };

  // ------------------------------------------------------------------
  // Convenience accessors for existing meson sectors
  // ------------------------------------------------------------------
  namespace PionFF {
    ResonanceParameters& Rho();
    ResonanceParameters& Rho1450();
    ResonanceParameters& Rho1700();
    ResonanceParameters& Rho2150();
    ResonanceParameters& Omega782();
    ResonanceParameters& Phi1020();
    void ReadAll();
  } // namespace PionFF

  namespace KaonFF {
    ResonanceParameters& Rho();
    ResonanceParameters& Rho1450();
    ResonanceParameters& Rho1700();
    ResonanceParameters& Rho2150();
    ResonanceParameters& Omega();
    ResonanceParameters& Omega1420();
    ResonanceParameters& Omega1650();
    ResonanceParameters& Phi();
    ResonanceParameters& Phi1680();
    void ReadAll();
  } // namespace KaonFF

  namespace ThreePionFF {
    ResonanceParameters& Omega782();
    ResonanceParameters& Phi1020();
    ResonanceParameters& Jpsi();
    void ReadAll();
  } // namespace ThreePionFF

  namespace TwoPhotonFF {
    // π0
    ResonanceParameters& Pi0Rho();
    ResonanceParameters& Pi0Omega();
    ResonanceParameters& Pi0Phi();
    // η
    ResonanceParameters& EtaRho();
    ResonanceParameters& EtaOmega();
    ResonanceParameters& EtaPhi();
    // η'
    ResonanceParameters& EtapRho();
    ResonanceParameters& EtapOmega();
    ResonanceParameters& EtapPhi();
    void ReadAll();
  } // namespace TwoPhotonFF

  // ------------------------------------------------------------------
  // Baryon sector – placeholder for future baryon form factors
  // (e.g., N* and Δ* resonances)
  // ------------------------------------------------------------------
  namespace BaryonFF {
    // Examples – you can add more as needed
    ResonanceParameters& N1440(); // Roper
    ResonanceParameters& N1520();
    ResonanceParameters& N1535();
    ResonanceParameters& Delta1232();
    ResonanceParameters& Delta1600();
    // ... add all resonances you need
    void ReadAll();
  } // namespace BaryonFF

  // ------------------------------------------------------------------
  // Utility functions (same as before, but now use registry)
  // ------------------------------------------------------------------
  void FFTableHead(const std::string& title);
  void FFTableFoot(const std::string& label, const Complex& f0);
  void CheckFFNormalisation(const std::string& tag, const Complex& f0,
                            bool requireunity = true);

  extern const int s_ffwidth; // table width

} // namespace METOOLS

#endif