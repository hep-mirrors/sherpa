#ifndef Remnant_Base_H
#define Remnant_Base_H

#include "Beam_Spectra_Handler.H"
#include "Particle_List.H"
#include "Blob.H"
#include "My_Vector.H"

namespace SHERPA {

  class Remnant_Base {
  public:

    enum TypeID { Hadron   =  1,
		  Electron =  2,
		  Photon   =  3 };

  protected:
    
    TypeID       m_type;
    unsigned int m_beam;

    ATOOLS::Particle_List m_parton[3];

    ATOOLS::Particle *p_last[2];
    Remnant_Base     *p_partner;

    ATOOLS::Blob  *p_beamblob;
    ATOOLS::Vec4D  m_pbeam;

    bool m_active;

    double Lambda2(double sp,double sp1,double sp2);

  public:

    // constructor
    Remnant_Base(const TypeID _m_type,const unsigned int _m_beam);

    // destructor
    virtual ~Remnant_Base();

    // setters
    void SetPartner(SHERPA::Remnant_Base *_p_partner);
    void ExtractParton(ATOOLS::Particle *_m_parton);

    // getters
    const TypeID       Type() const;
    const unsigned int Beam() const;

    ATOOLS::Particle *Last();
    ATOOLS::Blob     *BeamBlob();

    // member functions
    virtual void Clear();

    virtual bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist) = 0;
    virtual bool AdjustKinematics();

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
