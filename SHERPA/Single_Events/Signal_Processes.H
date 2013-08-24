#ifndef SHERPA_Single_Events_Signal_Processes_H
#define SHERPA_Single_Events_Signal_Processes_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "PDF/Remnant/Remnant_Base.H"

namespace METOOLS {
  class Amplitude2_Tensor;
}

namespace SHERPA {

  class Signal_Processes: public Event_Phase_Handler {
  private:

    Matrix_Element_Handler *p_mehandler;
    PDF::Remnant_Base      *p_remnants[2];
    METOOLS::Amplitude2_Tensor* p_atensor;

    double m_overweight;
    int m_setcolors, m_cmode;

    bool FillBlob(ATOOLS::Blob_List *const bl,ATOOLS::Blob *const b);

  public:

    Signal_Processes(Matrix_Element_Handler *);

    ~Signal_Processes();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp(const size_t & mode=0);
    void Finish(const std::string &);

    inline Matrix_Element_Handler *GetMEHandler() const { return p_mehandler; }

  };

}

#endif

