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

  /*!
    Everything Griffin needs to evaluate ONE registered process: the two
    external PDG codes and the identical-particle normalisation. Built once,
    by MakeProcess(), and owned by that process's Griffin_Born or
    Griffin_Virtual -- NOT shared class-static state. Griffin is 2->2 only,
    so the incoming pair is (initial_pdg, bar(initial_pdg)) and the outgoing
    pair is (final_pdg, bar(final_pdg)).

    This used to be three class-static members on Griffin_Interface
    (m_inital, m_final, m_norm), overwritten by RegisterProcess() every time
    a new process was constructed. With more than one Griffin-backed process
    registered in the same run (e.g. e+e- -> mu+mu- and Bhabha together),
    whichever process registered last silently won for ALL of them: every
    other process's per-point evaluation read the wrong flavours. Moving
    this into a per-process value, passed explicitly into every Evaluate*
    call, makes that class of bug impossible to reintroduce.
  */
  struct GriffinProcess {
    int    initial_pdg = 11;
    int    final_pdg   = 13;
    double norm        = 1.;
  };

  /*!
    Sherpa-facing interface to the Griffin electroweak loop library
    (https://arxiv.org/abs/2211.16272). Exactly one instance exists per run,
    created by Sherpa's ME_Generator_Base factory when "Griffin" is named in
    ME_GENERATORS; Initialize() fixes the electroweak input scheme once, from
    the runcard, and every registered process evaluates against that one,
    shared set of Standard Model inputs via Instance().

    What stays shared vs what became per-process, and why:
      - The SM inputs (m_griffin) and the requested order (m_order) describe
        the ELECTROWEAK SCHEME for the whole run: there is exactly one
        EW_SCHEME per Sherpa run, so these correctly stay here, as private
        members of the single instance, reached through Instance().
      - The two external flavours and any identical-particle normalisation
        (GriffinProcess, above) describe ONE process and differ between
        processes, so they are no longer members of this class at all.
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

    //! The one Griffin_Interface for this run. THROWs if called before
    //! Sherpa's ME_Generator factory has constructed it (i.e. before
    //! ME_GENERATORS: [..., Griffin, ...] has been processed).
    static Griffin_Interface &Instance();

    //! Build the per-process context for a soon-to-be-constructed
    //! Griffin_Born or Griffin_Virtual. Pure functions of their argument --
    //! no shared state is read or written, so registering one process
    //! cannot affect another's context.
    static GriffinProcess MakeProcess(const PHASIC::Process_Info &pi);
    static GriffinProcess MakeProcess(const PHASIC::External_ME_Args &args);

    //! (Higher order) - (LO), i.e. the loop contribution Griffin_Virtual
    //! multiplies onto the dressed Born, in Griffin's own normalisation.
    void EvaluateLoop(const ATOOLS::Vec4D_Vector &momenta,
                       const GriffinProcess &proc,
                       METOOLS::DivArrD &virt) const;
    void EvaluateBorn(const ATOOLS::Vec4D_Vector &momenta,
                       const GriffinProcess &proc,
                       double &born) const;

  private:
    void RegisterDefaults();

    //! s and cost = cos(theta) for a 2->2 point, in the convention every
    //! Evaluate* below shares: t = (p0-p2)^2, cost = 1 + 2t/s. Factored out
    //! so a future evaluator cannot independently reinvent (and get wrong)
    //! this kinematic mapping -- exactly what happened once already, when
    //! EvaluateBorn computed cost from an unrelated quantity. A named
    //! struct rather than std::pair<double,double>: it self-documents
    //! which field is which, and (unlike a structured binding) is capturable
    //! by value or reference in every lambda below without relying on the
    //! C++20 extension clang otherwise silently allows here.
    struct Kinematics { double s, cost; };
    static Kinematics KinVars(const ATOOLS::Vec4D_Vector &momenta);

    //! Differential-cross-section term shared by the LO, NLO and NNLO
    //! evaluators: build the vertex form factors and vector/axial-vector
    //! matrix element in whichever concrete FA/SW/Mat classes are given,
    //! and combine them into the standard helicity sum. Templated so the
    //! physics differs only in which of Griffin's own classes are passed in
    //! -- not in three independently-maintained copies of one formula.
    template <class FA, class SW, class Mat, class SM>
    static double DiffXSecTerm(int initial_pdg, int final_pdg,
                                double s, double cost, const SM &gx);

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
