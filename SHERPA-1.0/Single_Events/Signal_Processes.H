#ifndef Signal_Processes_h
#define Signal_Processes_h

#include "Event_Phase_Handler.H"
#include "Matrix_Element_Handler.H"
#include "Sample_Analysis.H"

namespace SHERPA {
  class Signal_Processes : public Event_Phase_Handler {
  private :
    Matrix_Element_Handler * p_mehandler;
    void FillBlob(APHYTOOLS::Blob *);
  public :
    Signal_Processes(Matrix_Element_Handler *);
    ~Signal_Processes();
    bool Treat(APHYTOOLS::Blob_List *, double &);
    void CleanUp();
  };
}
#endif

