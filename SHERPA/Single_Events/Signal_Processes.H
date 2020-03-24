#ifndef SHERPA_Single_Events_Signal_Processes_H
#define SHERPA_Single_Events_Signal_Processes_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "REMNANTS/Main/Remnant_Base.H"

namespace ATOOLS { class Scoped_Settings; }

namespace SHERPA {


  class Signal_Processes: public Event_Phase_Handler {
  private:

    Matrix_Element_Handler    *p_mehandler;
    REMNANTS::Remnant_Base    *p_remnants[2];

    double m_overweight;
    int m_cmode;
    bool m_setcolors, m_adddocumentation;

    bool FillBlob(ATOOLS::Blob_List* bl, ATOOLS::Blob* blob);
    int PrepareShowerBlob(ATOOLS::Blob_List* bl, ATOOLS::Blob* signalblob);

  public:

    Signal_Processes(Matrix_Element_Handler*);

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void CleanUp(const size_t & mode=0);
    void Finish(const std::string &);

    inline Matrix_Element_Handler *GetMEHandler() const { return p_mehandler; }

  };

}

#endif

