#ifndef Multiple_Interactions_H
#define Multiple_Interactions_H

#include "Event_Phase_Handler.H"
#include "MI_Handler.H"

namespace SHERPA {

  class Multiple_Interactions: public Event_Phase_Handler {
  private:

    double m_ecms, m_pperpmax, m_ecmsmax, m_xmax[2];

    MI_Handler *p_mihandler;

    bool CheckBlobList(const ATOOLS::Blob_List *bloblist);

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

