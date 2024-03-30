#ifndef REMNANTS_Main_Remnant_Handler_H
#define REMNANTS_Main_Remnant_Handler_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Kinematics_Generator.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include "REMNANTS/Tools/Beam_Decorrelator.H"
#include "PDF/Main/ISR_Handler.H"
#include "YFS/Main/YFS_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Return_Value.H"

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
    PDF::isr::id         m_id;
    Remnant_Base       * p_remnants[2];
    size_t               m_tags[2];
    strat::code          m_type;
    Kinematics_Generator m_kinematics;
    Colour_Generator     m_colours;
    Beam_Decorrelator    m_decorrelator;
    ATOOLS::Blob       * p_softblob;

    std::set<ATOOLS::Blob *> m_treatedshowerblobs;

    bool m_check, m_output;
    long int m_fails;
    
    void InitializeRemnants(PDF::ISR_Handler * isr,
          YFS::YFS_Handler *yfs,
			    BEAM::Beam_Spectra_Handler * beam,
			    const std::vector<size_t> & tags);
    void DefineRemnantStrategy();
    void InitializeKinematicsAndColours();
    bool CheckBeamBreakup(ATOOLS::Blob_List * bloblist);
    void InitBeamAndSoftBlobs(ATOOLS::Blob_List *const blobs,
			      const bool & isrescatter);
    ATOOLS::Blob_List::iterator FindInsertPositionForRescatter(ATOOLS::Blob_List *const bloblist,
							       const bool & isrescatter);
 public:
    Remnant_Handler(PDF::ISR_Handler * isr, YFS::YFS_Handler *yfs, BEAM::Beam_Spectra_Handler * beam,
		    const std::vector<size_t> & tags);
    ~Remnant_Handler();
    
    bool ExtractShowerInitiators(ATOOLS::Blob *const showerblob);
    bool Extract(ATOOLS::Particle * part,const unsigned int beam);
    void ConnectColours(ATOOLS::Blob *const showerblob);
    void Reset();
    
    ATOOLS::Return_Value::code MakeBeamBlobs(ATOOLS::Blob_List *const blobs,
					     ATOOLS::Particle_List *const=nullptr,
					     const bool & isrescatter=false);

    void SetImpactParameter(const double & b);
    
    inline bool NeedsToMakeBeamBlobs() {
      return (m_colours.HasColours(0) || m_colours.HasColours(0));
    }
    inline void SetScale2(const double & scale2) {
      for (size_t i=0;i<2;i++) p_remnants[i]->SetScale2(scale2);
    }
    inline const PDF::isr::id Id()  const { return m_id; }
    inline const strat::code Type() const { return m_type; }
    
    inline Remnant_Base * GetRemnant(const unsigned short int & beam) const {
      return (beam<2)?p_remnants[beam]:nullptr;
    }
    inline Primordial_KPerp * GetKPerp() { return m_kinematics.GetKPerp(); }
  };
  
  typedef std::map<PDF::isr::id, Remnant_Handler*> Remnant_Handler_Map;
}

#endif
