#ifndef Remnant_Base_H
#define Remnant_Base_H

#include "Remnant_Info.H"
#include "Particle_List.H"
#include "Blob.H"
#include <vector>

namespace SHERPA {

  class MI_Handler;

  struct rtp {
    
    enum code { 
      intact      =  0,
      qcd_remnant =  1,
      hadron      =  3,
      photon      =  5, 
      electron    =  8
    };

  };//end of struct rtp

  inline bool operator&(const rtp::code code1,const rtp::code code2)
  { return (rtp::code)((int)code1&(int)code2); }

  std::ostream &operator<<(std::ostream &ostr,const rtp::code code);

  class Remnant_Base: public ATOOLS::Remnant_Info {
  protected:
    
    rtp::code    m_type;
    unsigned int m_beam, m_errors;

    ATOOLS::Particle_List m_parton[2];

    ATOOLS::Particle *p_last[2];
    Remnant_Base     *p_partner;

    MI_Handler *p_mihandler;

    ATOOLS::Blob  *p_beamblob;
    ATOOLS::Vec4D  m_pbeam;

    double m_erem, m_pzrem, m_deltae, m_deltap;

    bool m_active, m_initialized;

    bool FindHardProcess(ATOOLS::Particle *const initiator,ATOOLS::Blob *&process);
    void FindIncoming(ATOOLS::Blob *const blob,const size_t beam,ATOOLS::Particle *&incoming);
    bool AcquireMass(const ATOOLS::Particle *left,
		     const ATOOLS::Particle *right,const double newsp);

    bool AdjustEnergy();

  public:

    // constructor
    Remnant_Base(const rtp::code type,const unsigned int beam);

    // destructor
    virtual ~Remnant_Base();

    // member functions
    virtual void Clear();

    virtual double MinimalEnergy(const ATOOLS::Flavour &flavour);

    virtual bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist) = 0;
    virtual bool AdjustKinematics();

    virtual void UnDo(); 

    // inline functions 
    inline void SetPartner(SHERPA::Remnant_Base *partner)   { p_partner=partner;     }
    inline void SetMIHandler(SHERPA::MI_Handler *mihandler) { p_mihandler=mihandler; }

    inline ATOOLS::Particle *Last()     { m_active=false; return p_last[0]; }
    inline ATOOLS::Blob     *BeamBlob() { return p_beamblob;                }

    inline rtp::code    Type() const { return m_type; }
    inline unsigned int Beam() const { return m_beam; }

    inline MI_Handler *const MIHandler() { return p_mihandler; }

    inline unsigned int Errors() const { return m_errors; }

    inline void ClearErrors() { m_errors=0; }

    inline void Extract(ATOOLS::Particle *parton) { m_parton[1].push_back(parton); }

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
