#ifndef Multiple_Interactions_H
#define Multiple_Interactions_H

#include "Event_Phase_Handler.H"
#include "MI_Handler.H"

namespace ATOOLS { class Jet_Finder;   }
namespace PDF    { class Remnant_Base; }

namespace SHERPA {

  class Multiple_Interactions: public Event_Phase_Handler {
  private:

    double m_ecms, m_ptmax, m_lastptmax, m_ecmsmax, m_emax[2];
    double m_weight, m_ntrials;

    MI_Handler         *p_mihandler;
    ATOOLS::Jet_Finder *p_jetfinder;
    PDF::Remnant_Base  *p_remnants[2];

    ATOOLS::Blob_List *p_bloblist;

    bool m_vetoed, m_diced;

    ATOOLS::Return_Value::code CheckBlobList(ATOOLS::Blob_List *const bloblist);
    bool                       VetoHardProcess(ATOOLS::Blob *const blob);
    
  public:

    // constructor
    Multiple_Interactions(MI_Handler *);

    // destructor
    ~Multiple_Interactions();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);

  }; // end of class Multiple_Interactions

} // end of namespace SHERPA

#endif

