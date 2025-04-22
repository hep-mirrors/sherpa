#ifndef Griffin_Interface_H
#define Griffin_Interface_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include <algorithm>
#include <sys/stat.h>

#include "classes.h"
#include "EWPOZ2.h"
#include "xscnnlo.h"
#include "xscenlo.h"
#include "SMval.h"


using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;
using namespace griffin;

namespace Griffin {
  struct griffinorder {
  enum code {
    lo       = 0,
    nlo      = 1,
    nnlo     = 2,
    nloe     =3,
  };
};

  std::istream &operator>>(std::istream &str, griffinorder::code &sm);
  std::ostream &operator<<(std::ostream &str, const griffinorder::code &sm);
  class Griffin_Interface: public PHASIC::ME_Generator_Base {
  public:
    Griffin_Interface();
    ~Griffin_Interface();
    bool Initialize(MODEL::Model_Base *const model,
          BEAM::Beam_Spectra_Handler *const beam,
          PDF::ISR_Handler *const isr,
          YFS::YFS_Handler *const yfs);
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    void RegisterDefaults();
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add)
    { return NULL; }
    int  PerformTests();
    bool NewLibraries() { return false; }
    void PrintLogo(std::ostream &s);
    static void RegisterProcess(const PHASIC::Process_Info &pi);
    static void RegisterProcess(const External_ME_Args& args);
    static void EvaluateLoop(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt);
    static void EvaluateLO(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt);
    static void EvaluateNLO(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt);
    static void EvaluateNLOE(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt);
    static void EvaluateNNLO(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt);
    static void EvaluateBorn(const Vec4D_Vector& momenta, double &born);
    static griffin::inval m_griffin;
    static griffinorder::code m_order;
    static int m_inital;
    static int m_final;
    static double m_norm;
    static double m_B;
  };

}
#endif