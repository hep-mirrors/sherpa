#ifndef PDF_Remnant_Remnant_Base_H
#define PDF_Remnant_Remnant_Base_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob.H"
#include <vector>

namespace BEAM { class Beam_Base; }

namespace PDF {

  struct rtp {
    
    enum code { 
      intact      =  0,
      qcd_remnant =  1,
      hadron      =  3,
      photon      =  5, 
      electron    =  8
    };

  };// end of struct rtp

  inline bool operator&(const rtp::code code1,const rtp::code code2)
  { return (rtp::code)((int)code1&(int)code2); }

  std::ostream &operator<<(std::ostream &ostr,const rtp::code code);

  class Remnant_Base {
  private:

    static std::set<ATOOLS::Particle*> s_last[2];

  protected:
    
    BEAM::Beam_Base *p_beam;

    rtp::code    m_type;
    unsigned int m_beam, m_errors;

    ATOOLS::Particle_List m_extracted, m_companions;

    ATOOLS::Particle *p_last[2];
    Remnant_Base     *p_partner;

    ATOOLS::Blob  *p_beamblob;
    ATOOLS::Vec4D  m_pbeam;

    double m_erem, m_pzrem, m_lastemin, m_deltae, m_deltap;
    double m_emin, m_scale, m_laste;

    bool m_active, m_initialized;

    bool FindLast(const short unsigned int side);
    bool ChooseLast();

  public:

    // constructor
    Remnant_Base(const rtp::code type,const unsigned int beam);

    // destructor
    virtual ~Remnant_Base();

    // member functions
    virtual void Clear();
    virtual void QuickClear();

    virtual double          MinimalEnergy(const ATOOLS::Flavour &flavour);

    virtual bool FillBlob(ATOOLS::Blob *beamblob,
			  ATOOLS::Particle_List *particlelist) = 0;
    virtual bool AdjustKinematics();
    virtual bool AdjustColors();

    virtual void UnDo(); 

    bool TestExtract(const ATOOLS::Flavour &flav,
		     const ATOOLS::Vec4D &mom) ;
    bool TestExtract(ATOOLS::Particle *parton);
    bool Extract(ATOOLS::Particle *parton);

    // inline functions 
    inline void SetPartner(PDF::Remnant_Base *partner) { p_partner=partner; }

    inline ATOOLS::Particle *Last()     { m_active=false; return p_last[0]; }
    inline ATOOLS::Blob     *BeamBlob() { return p_beamblob;                }

    inline rtp::code    Type() const { return m_type; }
    inline unsigned int Beam() const { return m_beam; }

    inline BEAM::Beam_Base *GetBeam() const       { return p_beam; }
    inline double           RemnantEnergy() const { return m_erem; }

    inline unsigned int Errors() const { return m_errors; }

    inline void ClearErrors() { m_errors=0; }

    inline void SetBeam(BEAM::Beam_Base *const beam) { p_beam=beam;   }
    inline void SetScale(const double scale)         { m_scale=scale; }

  };// end of class Remnant_Base

}// end of namespace PDF

#endif
