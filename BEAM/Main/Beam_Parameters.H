#ifndef BEAM_Main_Beam_Parameters_H
#define BEAM_Main_Beam_Parameters_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Info_Key.H"

namespace BEAM {
  enum class beammode {
    unknown              = 0,
    relic_density        = 1,
    collider             = 2,
    DM_annihilation      = 10,
    Fixed_Target         = 12
  };
  std::ostream& operator<<(std::ostream&,const beammode mode);

  enum class beamspectrum {
    monochromatic        = 0,
    Gaussian             = 1,
    laser_backscattering = 2,
    simple_Compton       = 3,
    EPA                  = 4,
    Pomeron              = 5,
    DM                   = 10,
    Fixed_Target         = 12,
    spectrum_reader      = 20,
    unknown              = 99
  };
  std::ostream& operator<<(std::ostream&,const beamspectrum spect);

  enum class collidermode {
    monochromatic = 0,
    spectral_1    = 1,
    spectral_2    = 2,
    both_spectral = 3,
    unknown       = 99
  };
  std::ostream& operator<<(std::ostream&,const collidermode mode);


  class Beam_Base;
  class Beam_Parameters {
  private:
    beammode m_beammode;
    std::array<beamspectrum, 2> m_beamspec;

    ATOOLS::Settings & m_settings;

    void RegisterDefaults();
    void RegisterDefaultBeams();
    void RegisterLaserDefaults();
    void RegisterEPADefaults();
    void RegisterPomeronDefaults();
    void RegisterDarkMatterDefaults();
    void InitializeFlav(kf_code flav);
    Beam_Base * InitializeSimpleCompton(int);
    /*!
      Characteristic parameters for laser back scattering are read in and the spectrum is
      instantiated through a call to
      Laser_Backscattering::Laser_Backscattering(const ATOOLS::Flavour,const double,const double,
                                                 const double,const double,const int,
                                                 const int,const int,bool &)
    */
    Beam_Base * InitializeLaserBackscattering(int);
    Beam_Base * InitializeSpectrumReader(int);
    /*!
      Characteristic parameters for monochromatic beams are read in and the spectrum is
      instantiated through a call to
      Monochromatic::Monochromatic(const ATOOLS::Flavour,const double,const double,bool &)
    */
    Beam_Base * InitializeMonochromatic(int);
    Beam_Base * InitializeEPA(int);
    Beam_Base * InitializePomeron(int);
    Beam_Base * InitializeDM_beam(int);
    Beam_Base * InitializeFixed_Target(int);
  public:
    Beam_Parameters();
    ~Beam_Parameters();

    bool SpecifyMode();
    /*!
      SpecifySpectra reads in the specific spectra of the beams and the corresponding generator
      and initializes the spectra through appropriate InitializeXXX methods. At the moment, two
      such methods are available, namely InitializeLaserBackscattering(int)
      and InitializeMonochromatic(int). For every new spectrum, such
      a method has to be provided as well and SpecifySpectra has to be modified accordingly.
      The idea is that in the long run, beam generators might be linked from outside by
      experimenters aiming at higher precision/better fit to beam data.
    */
    bool SpecifySpectra();
    Beam_Base * InitSpectrum(const size_t & num);

    const ATOOLS::Flavour GetFlavour(const std::string & tag,const size_t & pos=-1);
    const std::string     String(const std::string & tag,const int & pos=-1) const;
    const double          operator()(const std::string & tag,const int & pos=-1) const;
    const int             Switch(const std::string & tag,const int & pos=-1) const;
    const bool            On(const std::string & tag) const;

    inline beammode     GetMode() const { return m_beammode; }
    inline beamspectrum GetSpectrum(const size_t & pos) const {
      if (pos==0 || pos==1) return m_beamspec[pos];
      return beamspectrum::unknown;
    }
  };
}

#endif

