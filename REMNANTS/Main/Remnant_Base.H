#ifndef REMNANTS_Main_Remnant_Base_H
#define REMNANTS_Main_Remnant_Base_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob.H"
#include <vector>

namespace BEAM { class Beam_Base; }

namespace REMNANTS {  
  struct rtp {
    enum code {
      none   = 0,
      intact = 1,
      lepton = 2,
      hadron = 3,
      photon = 4 
    };
  };

  inline bool operator&(const rtp::code code1,const rtp::code code2) {
    return (rtp::code)((int)code1&(int)code2);
  }
  std::ostream &operator<<(std::ostream &ostr,rtp::code code);
}
  
namespace REMNANTS {
  typedef std::list<ATOOLS::Flavour> FlavourList;

  class Colour_Generator;
  
  class Remnant_Base {
  protected:
    rtp::code    m_type;
    unsigned int m_beam;

    BEAM::Beam_Base  * p_beam;
    ATOOLS::Blob     * p_beamblob;
    ATOOLS::Part_List  m_extracted, m_spectators;
    Colour_Generator * p_colours;

    double m_energy, m_scale2;

    long unsigned int m_errors;
  public:
    Remnant_Base(rtp::code type,unsigned int beam);
    virtual ~Remnant_Base();

    virtual bool Extract(ATOOLS::Particle *parton);
    virtual bool TestExtract(ATOOLS::Particle *parton);
    virtual bool TestExtract(const ATOOLS::Flavour &flav,
			     const ATOOLS::Vec4D &mom)=0;
    virtual ATOOLS::Blob * MakeBlob();
    virtual bool FillBlob(ATOOLS::ParticleMomMap *ktmap=nullptr,
			  const bool & copy=true)=0;
    virtual void MakeSpectator(ATOOLS::Particle * parton) {}
    
    virtual void Reset(const bool & DIS=false);
    inline rtp::code    Type() const                 { return m_type; }
    inline unsigned int Beam() const                 { return m_beam; }
    inline void SetBeam(BEAM::Beam_Base * beam)      { p_beam = beam; }
    inline BEAM::Beam_Base * GetBeam() const         { return p_beam; }
    inline ATOOLS::Blob * GetBlob()                  { return p_beamblob; }
    inline const ATOOLS::Vec4D & InMomentum() const  { return p_beam->OutMomentum(); }
    inline void SetScale2(const double & scale2)     { m_scale2 = scale2; }
    inline ATOOLS::Part_List * GetExtracted()        { return &m_extracted; }
    inline ATOOLS::Part_List * GetSpectators()       { return &m_spectators; }
    inline virtual ATOOLS::Particle * GetRecoiler()  { return nullptr; }
    inline virtual ATOOLS::Particle * GetSpectator() { return nullptr; }
    inline void SetColours(Colour_Generator * cgen)  { p_colours = cgen; } 
  };
}

#endif
