#ifndef SHERPA_Single_Events_Multiple_Interactions_H
#define SHERPA_Single_Events_Multiple_Interactions_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"

namespace SHERPA {

  class Multiple_Interactions: public Event_Phase_Handler {
  private:
    double m_hardveto, m_ptmax_fac;
    double m_ecms, m_ptmax, m_lastptmax, m_ecmsmax, m_emax[2];
    double m_defaultptmax;
    double m_weight, m_ntrials;
    bool   m_isfirstMB;

    MI_Handler_Map    * p_mihandlers;
    MI_Handler        * p_activeMI;
    ATOOLS::Blob_List * p_bloblist;
    ATOOLS::Blob      * p_lastblob;

    std::set<ATOOLS::Blob *> m_ISblobs;
    bool m_vetoed, m_newevent;

    ATOOLS::Return_Value::code m_result;
    void Reset();
    void SwapToRescatter();
    void ResetMIHandler();
    bool CheckBlobList();    
    bool InitNewEvent();
    void InitMinBiasEvent();
    bool BeamsViable();
    bool ExtractISInfo(ATOOLS::Blob * blob);
    void SwitchPerturbativeInputsToMIs();
    bool MIKinematics();

  public:
    Multiple_Interactions(MI_Handler_Map *);
    ~Multiple_Interactions();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void CleanUp(const size_t & mode=0);
    void Finish(const std::string &);
  };
}

#endif

