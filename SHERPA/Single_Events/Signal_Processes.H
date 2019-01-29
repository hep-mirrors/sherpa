#ifndef SHERPA_Single_Events_Signal_Processes_H
#define SHERPA_Single_Events_Signal_Processes_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Phys/Variations.H"

namespace ATOOLS { class Scoped_Settings; }

namespace SHERPA {


  class Signal_Processes: public Event_Phase_Handler {
  private:

    Matrix_Element_Handler    *p_mehandler;
    ATOOLS::Variations        *p_variations;
    ATOOLS::Variation_Weights  m_variationweights;
    REMNANTS::Remnant_Base    *p_remnants[2];

    double m_overweight;
    int m_cmode;
    bool m_setcolors, m_adddocumentation;

    bool FillBlob(ATOOLS::Blob_List *const bl,ATOOLS::Blob *const b);

  public:

    Signal_Processes(Matrix_Element_Handler*,
                     ATOOLS::Variations*);

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp(const size_t & mode=0);
    void Finish(const std::string &);

    inline Matrix_Element_Handler *GetMEHandler() const { return p_mehandler; }

  };

}

#endif

