#ifndef ATOOLS_Phys_Fastjet_Helpers_H
#define ATOOLS_Phys_Fastjet_Helpers_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Phys/fjcore.hh"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

#include <vector>
#include <limits>

namespace ATOOLS {
  struct Jet_Input;
  class  Jet_Identification;
  typedef std::vector<Jet_Input> Jet_Inputs;
  typedef std::vector<Jet_Identification *> Jet_Identifications;

  bool ToBeClustered(const ATOOLS::Flavour& flav, int bmode);
  bool ToBeClustered(const ATOOLS::Flavour& flav,
                     const ATOOLS::Vec4D& mom,
                     const ATOOLS::Jet_Inputs& jetinputs);
  fjcore::PseudoJet MakePseudoJet(const Flavour& flav, const Vec4D& mom);

  bool BTag(const fjcore::PseudoJet& jet, int bmode);
  Flavour FlavourTag(const fjcore::PseudoJet& jet,
                     const Jet_Identifications& jetids,
                     const kf_code& notagkf);


  struct Jet_Input {
  public:
    ATOOLS::Flavour m_fl;
    double          m_ptmin,m_etamax,m_ymax;
    Jet_Input(const ATOOLS::Flavour& fl=ATOOLS::Flavour(kf_none),
              const double& ptmin=0.,
              const double& etamax=std::numeric_limits<double>::max(),
              const double& ymax=std::numeric_limits<double>::max()) :
      m_fl(fl), m_ptmin(ptmin), m_etamax(etamax), m_ymax(ymax) {}
  };

  std::ostream &operator<<(std::ostream &str,const Jet_Input &ji);

  struct JetIdMode {
    enum code {
      unknown  = 0,
      larger   = 1,
      relative = 2
    };
  };

  inline JetIdMode::code operator|(const JetIdMode::code idm1,
                                   const JetIdMode::code idm2)
  { return (JetIdMode::code)((int)idm1|(int)idm2); }
  inline const JetIdMode::code &operator|=(JetIdMode::code &idm1,
                                           const JetIdMode::code idm2)
  { return idm1=(JetIdMode::code)((int)idm1|(int)idm2); }
  inline JetIdMode::code operator&(const JetIdMode::code idm1,
                                   const JetIdMode::code idm2)
  { return (JetIdMode::code)((int)idm1&(int)idm2); }
  inline const JetIdMode::code &operator&=(JetIdMode::code &idm1,
                                           const JetIdMode::code idm2)
  { return idm1=(JetIdMode::code)((int)idm1&(int)idm2); }

  std::ostream &operator<<(std::ostream &str,const JetIdMode::code &idm);


  class Jet_Identification {
  private:
    ATOOLS::Flavour  m_flid;
    double           m_ptmin, m_etmin, m_emin;
    JetIdMode::code  m_jetidmode;
  public:
    Jet_Identification(const ATOOLS::Flavour& flid, const double& ptmin,
                       const double& etmin, const double& emin,
                       const JetIdMode::code& jetidmode);
    ~Jet_Identification() {}

    inline const ATOOLS::Flavour& Flavour()   const { return m_flid; }
    inline const double& PTMin()              const { return m_ptmin; }
    inline const double& ETMin()              const { return m_etmin; }
    inline const double& EMin()               const { return m_emin; }
    inline const JetIdMode::code& JetIDMode() const { return m_jetidmode; }
  };

  std::ostream &operator<<(std::ostream &str,const Jet_Identification &jid);

}

#endif
