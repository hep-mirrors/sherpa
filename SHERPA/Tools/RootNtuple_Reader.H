#ifndef SHERPA_Tools_RootNtuple_Reader_H
#define SHERPA_Tools_RootNtuple_Reader_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "PHASIC++/Process/NLO_Helpers.H"
#include <vector>
#ifdef USING__ROOT
#include "TChain.h"
#endif

namespace SHERPA {
  class RootNtuple_Reader : public Event_Reader_Base {
#ifdef USING__ROOT
    static const Int_t s_kMaxParticle = 10;
    UInt_t m_id;
    Int_t m_nparticle;
    Float_t p_px[s_kMaxParticle];
    Float_t p_py[s_kMaxParticle];
    Float_t p_pz[s_kMaxParticle];
    Float_t p_E[s_kMaxParticle];
    Int_t p_kf[s_kMaxParticle];

    Double_t m_wgt,m_wgt2;
    TChain* p_f;
#endif

    int                m_eventmode;
    size_t             m_evtid,m_subevtid,m_evtcnt,m_entries,m_evtpos;
    PHASIC::NLO_subevtlist m_nlos;

    void               CloseFile();
    bool               ReadInEntry();
    bool               ReadInFullEvent(ATOOLS::Blob_List *);
    bool               ReadInSubEvent(ATOOLS::Blob_List *);
  public:
    RootNtuple_Reader(const std::string & =std::string(""),const std::string & =std::string(""),int mode=0);
    bool        FillBlobs(ATOOLS::Blob_List *);
  };
}

#endif
