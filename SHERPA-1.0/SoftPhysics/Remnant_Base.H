#ifndef Remnant_Base_H
#define Remnant_Base_H

#include "Remnant_Info.H"
#include "Particle_List.H"
#include "Blob.H"
#include <vector>

namespace SHERPA {

  class Remnant_Base: public ATOOLS::Remnant_Info {
  public:

    enum TypeID { 
      Intact   =  0,
      Hadron   =  1,
      Electron =  2,
      Photon   =  3 
    };

  protected:
    
    TypeID       m_type;
    unsigned int m_beam, m_errors;

    ATOOLS::Particle_List m_parton[3];

    ATOOLS::Particle *p_last[2];
    Remnant_Base     *p_partner;

    ATOOLS::Blob  *p_beamblob;
    ATOOLS::Vec4D  m_pbeam;

    double m_erem, m_pzrem, m_deltae, m_deltap;

    bool m_active, m_initialized;

    void FindHardProcess(ATOOLS::Particle *const initiator,ATOOLS::Blob *&process);
    void FindIncoming(ATOOLS::Blob *const blob,const size_t beam,ATOOLS::Particle *&incoming);
    bool AcquireMass(ATOOLS::Blob *const process,const double newsp);

    bool AdjustEnergy();

  public:

    // constructor
    Remnant_Base(const TypeID type,const unsigned int beam);

    // destructor
    virtual ~Remnant_Base();

    // member functions
    void Clear();
    void DeleteRemnants();

    virtual double MinimalEnergy(const ATOOLS::Flavour &flavour);

    virtual bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist) = 0;
    virtual bool AdjustKinematics();

    virtual void UnDo(); 

    // inline functions 
    inline void SetPartner(SHERPA::Remnant_Base *_p_partner) { p_partner=_p_partner; }

    inline ATOOLS::Particle *Last()     { m_active=false; return p_last[0]; }
    inline ATOOLS::Blob     *BeamBlob() { return p_beamblob;                }

    inline TypeID       Type() const { return m_type; }
    inline unsigned int Beam() const { return m_beam; }

    inline unsigned int Errors() const { return m_errors; }

    inline void ClearErrors() { m_errors=0; }

    inline void ExtractParton(ATOOLS::Particle *parton) 
    { if (parton!=NULL) m_parton[1].push_back(parton); }

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
