#ifndef SHERPA_Tools_RootNtuple_Reader_H
#define SHERPA_Tools_RootNtuple_Reader_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include <vector>
#include <string.h>

namespace PDF { class ISR_Handler; }

namespace MODEL { class Model_Base; }

namespace PHASIC { class Scale_Setter_Base; }

namespace SHERPA {

  struct RootNTupleReader_Variables;

  struct RR_Process_Info {
    char m_type[2];
    ATOOLS::Flavour_Vector m_fl;
    inline RR_Process_Info(char *type,int n,ATOOLS::Flavour *fl):
      m_fl(ATOOLS::Flavour_Vector(fl,&fl[n]))
    { strcpy(m_type,type); }
    bool operator<(const RR_Process_Info pi) const;
  };

  std::ostream &operator<<(std::ostream &s,const RR_Process_Info &pi);

  class RootNtuple_Reader : public Event_Reader_Base {
    RootNTupleReader_Variables *p_vars;
    int                m_check,m_calc,m_ftype,m_lomode;
    size_t             m_evtid,m_subevtid,m_evtcnt,m_entries,m_evtpos;
    ATOOLS::NLO_subevtlist m_nlos;
    ATOOLS::PDF_Info   m_pdfinfo;
    ATOOLS::ME_Weight_Info m_mewgtinfo;
    PHASIC::Scale_Setter_Arguments m_sargs;
    PHASIC::KFactor_Setter_Arguments m_kargs;
    std::map<RR_Process_Info,PHASIC::Process_Base*> m_procs;
    PDF::ISR_Handler *p_isr;
    double m_ecms, m_xf1, m_xf2;
    std::string m_treename;

    void               CloseFile();
    bool               ReadInEntry();
    bool               ReadInFullEvent(ATOOLS::Blob_List *);
    double CalculateWeight(const double &mur2,const double &muf2,const int mode);
  public:
    RootNtuple_Reader(const Input_Arguments &args,int exact=0,int ftype=0);
    ~RootNtuple_Reader();
    bool        FillBlobs(ATOOLS::Blob_List *);
  };
  class ERootNtuple_Reader {};
  class EDRootNtuple_Reader {};
}

#endif
