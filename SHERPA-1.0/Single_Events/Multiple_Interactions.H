#ifndef Multiple_Interactions_H
#define Multiple_Interactions_H

#include "Event_Phase_Handler.H"
#include "Sample_Analysis.H"
#include "MI_Handler.H"

namespace SHERPA {

  class Multiple_Interactions : public Event_Phase_Handler {
  private :

    bool m_diced;
    double m_pperpmax, m_ecmsmax, m_xmax[2];

    MI_Handler *p_mihandler;

  public :

    // constructor
    Multiple_Interactions(MI_Handler *_p_mihandler);

    // destructor
    ~Multiple_Interactions();

    // member functions
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);

  }; // end of class Multiple_Interactions

} // end of namespace SHERPA

#endif

