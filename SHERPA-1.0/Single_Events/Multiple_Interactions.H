#ifndef Multiple_Interactions_H
#define Multiple_Interactions_H

#include "Event_Phase_Handler.H"
#include "MI_Handler.H"

namespace SHERPA {

  class Remnant_Base;

  class Multiple_Interactions: public Event_Phase_Handler {
  private:

    double m_ecms, m_pperpmax, m_ecmsmax, m_xmax[2];

    MI_Handler   *p_mihandler;
    Remnant_Base *p_remnants[2];

    ATOOLS::Blob_List *p_bloblist;

    std::set<ATOOLS::Blob*> m_deleted;

    bool m_vetoed;

    bool CheckBlobList(ATOOLS::Blob_List *const bloblist);
    bool VetoHardProcess(ATOOLS::Blob *const blob);
    
    void DeleteConnected(ATOOLS::Blob *const blob);

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

