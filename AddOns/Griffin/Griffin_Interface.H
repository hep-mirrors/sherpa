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
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include <variant>

#include "classes.h"
#include "EWPOZ2.h"
#include "xscnnlo.h"
#include "xscenlo.h"
#include "SMvalGMwMz.h"

// Deliberately no "using namespace" here: this header is included by three
// translation units (Griffin_Interface.C, Griffin_Virtual.C, Griffin_Born.C)
// plus anything that ever includes one of those, and a using-directive at
// header scope would leak PHASIC/MODEL/ATOOLS/std/griffin names into all of
// them. Each .C file brings in what it needs at file scope instead.

namespace Griffin {

  struct griffinorder {
    enum code {
      lo   = 0,
      nlo  = 1,
      nnlo = 2,
      nloe = 3,
    };
  };

  std::istream &operator>>(std::istream &str, griffinorder::code &sm);
  std::ostream &operator<<(std::ostream &str, const griffinorder::code &sm);


  struct GriffinProcess {
    int    initial_pdg = 11;
    int    final_pdg   = 13;
    double norm        = 1.;
  };

  /*!
    Sherpa-facing interface to the Griffin electroweak loop library
    (https://arxiv.org/abs/2211.16272). 
  */
  class Griffin_Interface: public PHASIC::ME_Generator_Base {
  public:
    Griffin_Interface();
    ~Griffin_Interface() override;

    bool Initialize(MODEL::Model_Base *const model,
                     BEAM::Beam_Spectra_Handler *const beam,
                     PDF::ISR_Handler *const isr,
                     YFS::YFS_Handler *const yfs) override;

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &,
                                             bool) override
    { return nullptr; }
    int  PerformTests() override { return 1; }
    bool NewLibraries() override { return false; }
    void PrintLogo(std::ostream &s);

    static Griffin_Interface &Instance();

    static GriffinProcess MakeProcess(const PHASIC::Process_Info &pi);
    static GriffinProcess MakeProcess(const PHASIC::External_ME_Args &args);

    //! (Higher order) - (LO), i.e. the loop contribution Griffin_Virtual
    //! multiplies onto the dressed Born, in Griffin's own normalisation.
    void EvaluateLoop(const ATOOLS::Vec4D_Vector &momenta,
                       const GriffinProcess &proc,
                       METOOLS::DivArrD &virt) const;
    bool HigherOrderHelicity(const ATOOLS::Vec4D_Vector &momenta,
                             const GriffinProcess &proc,
                             std::vector<Complex> &a0, std::vector<Complex> &a1,
                             std::vector<std::vector<int> > &hel) const;
    void EvaluateHigherOrder(const ATOOLS::Vec4D_Vector &momenta,
                             const GriffinProcess &proc,
                             METOOLS::DivArrD &virt) const;
    void EvaluateBorn(const ATOOLS::Vec4D_Vector &momenta,
                       const GriffinProcess &proc,
                       double &born) const;

  private:
    void RegisterDefaults();

    struct Kinematics { double s, cost; };
    static Kinematics KinVars(const ATOOLS::Vec4D_Vector &momenta);

    template <class FA, class SW, class Mat, class SM>
    static double DiffXSecTerm(int initial_pdg, int final_pdg,
                                double s, double cost, const SM &gx);

    //! the four vector/axial combinations (vv, av, va, aa) of one order
    template <class FA, class SW, class Mat, class SM>
    static void FormFactors(int initial_pdg, int final_pdg, double s, double cost,
                            const SM &gx, Complex res[4]);
    //! the axial sign conventions of HigherOrderHelicity (0: not yet fixed)
    mutable int m_axsign_i = 0, m_axsign_f = 0;
    double EvaluateLO  (const ATOOLS::Vec4D_Vector &momenta, const GriffinProcess &proc) const;
    double EvaluateNLO (const ATOOLS::Vec4D_Vector &momenta, const GriffinProcess &proc) const;
    double EvaluateNLOE(const ATOOLS::Vec4D_Vector &momenta, const GriffinProcess &proc) const;
    double EvaluateNNLO(const ATOOLS::Vec4D_Vector &momenta, const GriffinProcess &proc) const;

    static Griffin_Interface *p_instance;

    griffinorder::code      m_order    = griffinorder::nnlo;
    MODEL::ew_scheme::code  m_ewscheme = MODEL::ew_scheme::alphamZ;
    std::variant<griffin::SMval, griffin::SMvalGMwMz> m_griffin;
  };

}
#endif
