#ifndef REMNANTS_Main_Remnant_Base_H
#define REMNANTS_Main_Remnant_Base_H

#include "REMNANTS/Tools/Form_Factor.H"
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
    ATOOLS::Flavour    m_beamflav;
    rtp::code          m_type;
    size_t             m_beam, m_tag;
    BEAM::Beam_Base  * p_beam;
    Form_Factor      * p_ff;
    ATOOLS::Blob     * p_beamblob;
    ATOOLS::Vec4D      m_position;
    ATOOLS::Part_List  m_extracted, m_spectators;
    Colour_Generator * p_colours;

    double m_residualE, m_scale2;
    
    rtp::code FixType(ATOOLS::Flavour & flav);
  public:
    Remnant_Base(const ATOOLS::Flavour & flav,const size_t & beam,
		 const size_t & tag=0);
    virtual ~Remnant_Base();

    virtual bool Extract(ATOOLS::Particle *parton);
    virtual bool TestExtract(ATOOLS::Particle *parton);
    virtual bool TestExtract(const ATOOLS::Flavour &flav,
			     const ATOOLS::Vec4D &mom)=0;
    virtual ATOOLS::Blob * MakeBlob();
    virtual bool FillBlob(ATOOLS::ParticleMomMap *ktmap=nullptr,
			  const bool & copy=true)=0;
    virtual void MakeSpectator(ATOOLS::Particle * parton) {}
    virtual ATOOLS::Vec4D IncomingMomentum();
    
    virtual void Reset(const bool & resc=false,const bool & DIS=false);
    virtual void SetResidualEnergy(const double & E=-1.) {
      m_residualE = (E>=0.) ? E : p_beam->OutMomentum(m_tag)[0];
    }
    
    inline rtp::code                  Type() const       { return m_type; }
    inline size_t                     Beam() const       { return m_beam; }
    inline const ATOOLS::Vec4D      & Position() const   { return m_position; }
    inline BEAM::Beam_Base          * GetBeam() const    { return p_beam; }
    inline Form_Factor              * GetFF() const      { return p_ff; }
    inline ATOOLS::Blob             * GetBlob()          { return p_beamblob; }
    inline const ATOOLS::Vec4D      & InMomentum() const { return p_beam->OutMomentum(m_tag); }
    inline ATOOLS::Part_List        * GetExtracted()     { return &m_extracted; }
    inline ATOOLS::Part_List        * GetSpectators()    { return &m_spectators; }
    inline virtual ATOOLS::Particle * GetRecoiler()      { return nullptr; }
    inline virtual ATOOLS::Particle * GetSpectator()     { return nullptr; }
    inline void SetBeam(BEAM::Beam_Base * beam)          { p_beam = beam; }
    inline void SetScale2(const double & scale2)         { m_scale2 = scale2; }
    inline void SetColours(Colour_Generator * cgen)      { p_colours = cgen; }
    inline void SetPosition(const ATOOLS::Vec4D & pos)   { m_position = pos; }
  };
}

#endif
