#ifndef BEAM_Main_Beam_Parameters_H
#define BEAM_Main_Beam_Parameters_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Info_Key.H"

namespace BEAM {
  struct beammode {
    enum code {
      unknown              = 0,
      relic_density        = 1,
      collider             = 2,
      DM_annihilation      = 10
    };
  };
  std::ostream& operator<<(std::ostream&,const beammode::code mode);

  struct beamspectrum {
    enum code {
      monochromatic        = 0,
      Gaussian             = 1,
      laser_backscattering = 2,
      simple_Compton       = 3,
      EPA                  = 4,
      DM                   = 10,
      spectrum_reader      = 20,
      unknown              = 99
    };
  };
  std::ostream& operator<<(std::ostream&,const beamspectrum::code spect);

  struct collidermode {
    enum code {
      monochromatic = 0,
      spectral_1    = 1,
      spectral_2    = 2,
      both_spectral = 3,
      unknown       = 99
    };
  };
  std::ostream& operator<<(std::ostream&,const collidermode::code mode);

  
  class Beam_Base;
  class Beam_Parameters {
  private:
    beammode::code     m_beammode;
    beamspectrum::code m_beamspec[2];

    ATOOLS::Settings & m_settings;

    void RegisterDefaults();
    void RegisterDefaultBeams();
    void RegisterLaserDefaults();
    void RegisterEPADefaults();
    void RegisterDarkMatterDefaults();
    //void RegisterSpectrumReaderDefaults();
    void InitializeFlav(kf_code flav);
    void SetBeamMode();
    Beam_Base * InitializeSimpleCompton(int);
    Beam_Base * InitializeLaserBackscattering(int);
    Beam_Base * InitializeSpectrumReader(int);
    Beam_Base * InitializeMonochromatic(int);
    Beam_Base * InitializeEPA(int);
    Beam_Base * InitializeDM_beam(int);
  public:
    Beam_Parameters();
    ~Beam_Parameters();

    bool SpecifyMode();
    bool SpecifySpectra();
    Beam_Base * InitSpectrum(const size_t & num);

    const ATOOLS::Flavour GetFlavour(const std::string & tag,const size_t & pos=-1);
    const std::string     String(const std::string & tag,const int & pos=-1) const;
    const double          operator()(const std::string & tag,const int & pos=-1) const;
    const int             Switch(const std::string & tag,const int & pos=-1) const;
    const bool            On(const std::string & tag) const;
    
    const inline beammode::code     GetMode() const { return m_beammode; }
    const inline beamspectrum::code GetSpectrum(const size_t & pos) const {
      if (pos==0 || pos==1) return m_beamspec[pos];
      return beamspectrum::unknown; 
    }
  };
}

#endif

