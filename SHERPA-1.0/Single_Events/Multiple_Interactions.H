#ifndef Multiple_Interactions_H
#define Multiple_Interactions_H

#include "Event_Phase_Handler.H"
#include "MI_Handler.H"

namespace SHERPA {

  class Remnant_Base;
  class Matrix_Element_Handler;

  class Multiple_Interactions: public Event_Phase_Handler {
  private:

    double m_ecms, m_ptmax, m_lastptmax, m_ecmsmax, m_emax[2];
    double m_weight, m_ntrials;

    MI_Handler   *p_mihandler;
    Remnant_Base *p_remnants[2];

    Matrix_Element_Handler *p_mehandler;

    ATOOLS::Blob_List *p_bloblist;

    bool m_vetoed, m_diced;

    bool CheckBlobList(ATOOLS::Blob_List *const bloblist);
    bool VetoHardProcess(ATOOLS::Blob *const blob);
    
  public:

    // constructor
    Multiple_Interactions(MI_Handler *mihandler);

    // destructor
    ~Multiple_Interactions();

    // member functions
    bool Treat(ATOOLS::Blob_List *bloblist, double &weight);
    void CleanUp();
    void Finish(const std::string &resultpath);

  }; // end of class Multiple_Interactions

} // end of namespace SHERPA

#endif

