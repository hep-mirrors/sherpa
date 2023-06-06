#ifndef REMNANTS_Main_Remnant_Handler_H
#define REMNANTS_Main_Remnant_Handler_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Kinematics_Generator.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include "REMNANTS/Tools/Beam_Decorrelator.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Return_Value.H"

namespace PDF {
  class ISR_Handler;
}

namespace BEAM {
  class Beam_Spectra_Handler;
}

namespace REMNANTS {
  struct strat {
    enum code {
      none   = 0,
      simple = 1,
      ll     = 2,
      DIS1   = 4,
      DIS2   = 5,
      hh     = 8
    };
  };
  class Primordial_KPerp;

  class Remnant_Handler {
  private:
    Remnant_Base       * p_remnants[2];
    strat::code          m_type;
    Kinematics_Generator m_kinematics;
    Colour_Generator     m_colours;
    Beam_Decorrelator    m_decorrelator;
    ATOOLS::Blob       * p_softblob;

    std::set<ATOOLS::Blob *> m_treatedshowerblobs;

    bool m_active, m_check, m_output;
    
    void InitializeRemnants(PDF::ISR_Handler * isr,
			    BEAM::Beam_Spectra_Handler * beam);
    void DefineRemnantStrategy();
    void InitializeKinematicsAndColours();
    bool CheckBeamBreakup(ATOOLS::Blob_List * bloblist);
    void InitBeamAndSoftBlobs(ATOOLS::Blob_List *const blobs);
 public:
    Remnant_Handler(PDF::ISR_Handler * isr,BEAM::Beam_Spectra_Handler * beam);
    ~Remnant_Handler();
    
    bool ExtractShowerInitiators(ATOOLS::Blob *const showerblob);
    bool Extract(ATOOLS::Particle * part,const unsigned int beam);
    bool TestExtract(ATOOLS::Particle *parton,const unsigned int beam) {
      return !(m_active) || p_remnants[beam]->TestExtract(parton);
    }
    bool TestExtract(const ATOOLS::Flavour &flav,const ATOOLS::Vec4D &mom,const unsigned int beam) {
      return !(m_active) || p_remnants[beam]->TestExtract(flav,mom);
    }
    void ConnectColours(ATOOLS::Blob *const showerblob);
    void Reset();
    void Reset(const unsigned int beam) { p_remnants[beam]->Reset(); }
    ATOOLS::Return_Value::code MakeBeamBlobs(ATOOLS::Blob_List *const blobs,
					     ATOOLS::Particle_List *const=nullptr);

    void SetImpactParameter(const double & b);
    inline const strat::code Type()      const { return m_type; }
    inline const bool IsActive()         const { return m_active; }
    inline void SetActive(const bool & active) { m_active = active; } 
    inline void SetScale2(const double & scale2) {
      for (size_t i=0;i<2;i++) p_remnants[i]->SetScale2(scale2);
    }
    inline Remnant_Base * GetRemnant(const unsigned short int & beam) const {
      return (beam<2)?p_remnants[beam]:nullptr;
    }
    inline Primordial_KPerp * GetKPerp() { return m_kinematics.GetKPerp(); }
  };
}

#endif
