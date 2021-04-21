#ifndef SHERPA_SoftPhysics_Fragmentation_Handler_H
#define SHERPA_SoftPhysics_Fragmentation_Handler_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif

namespace AHADIC {
  class Ahadic;
}

namespace SHERPA {

  class Fragmentation_Handler: public ATOOLS::Terminator_Object {
  private:
    std::string m_fragmentationmodel;
    int         m_mode;
    bool        m_shrink, m_flagpartonics;

    AHADIC::Ahadic * p_ahadic;
#ifdef USING__PYTHIA
    Lund_Interface * p_lund;
#endif

    void RegisterDefaults();
    void Shrink(ATOOLS::Blob_List *);
  public:
    Fragmentation_Handler(std::string);
    ~Fragmentation_Handler();

    ATOOLS::Return_Value::code 
    operator()(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
    const std::string & FragmentationModel() const { return m_fragmentationmodel; }
    const int & Mode() const { return m_mode; }
#ifdef USING__PYTHIA
    inline Lund_Interface * GetLundInterface() const { return p_lund; }
#endif
    inline AHADIC::Ahadic * GetAhadic()        const { return p_ahadic; }
  };
}
#endif

