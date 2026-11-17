#include <catch2/catch_all.hpp>
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Phys/Selector_List.H"
#include <cmath>
#include <memory>
#include <string>

using namespace ATOOLS;
using namespace PHASIC;

namespace {

  void InitGlobalsOnce()
  {
    static bool initialised(false);
    if (initialised) return;
    // as in the Sherpa constructor, but without ATOOLS::exh, which would
    // replace the signal handlers of Catch2
    mpi = new My_MPI();
    msg = new Message();
    rpa = new Run_Parameter();
    // Decaydata.yaml is only in the share path after installation
    rpa->gen.SetVariable("SHERPA_SHARE_PATH",
                         std::string(SHERPA_SOURCE_DIR)+"/HADRONS++");
    static char name[] = "sherpa_tests";
    static char *argv[] = {name, NULL};
    Settings::InitializeMainSettings(1,argv);
    AddParticle(kf_gluon,.0,.0,.0,0,8,2,-1,1,1,0,"G","G","G","G");
    AddParticle(kf_photon,.0,.0,.0,0,0,2,-1,1,1,0,"P","P","\\gamma","\\gamma");
    initialised=true;
  }

  class Flavours_Only_Process: public Process_Base {
  public:
    Flavours_Only_Process(const Flavour_Vector &flavs,const size_t &nin)
    {
      m_flavs=flavs;
      m_nin=nin;
      m_nout=flavs.size()-nin;
    }
    void SetScale(const Scale_Setter_Arguments &) override {}
    void SetKFactor(const KFactor_Setter_Arguments &) override {}
    void InitializeTheReweighting(Variations_Mode) override {}
    size_t Size() const override { return 1; }
    Process_Base *operator[](const size_t &) override { return this; }
    Weight_Info *OneEvent(const int,Variations_Mode,const int) override
    { return NULL; }
    Weights_Map Differential(const Vec4D_Vector &,Variations_Mode) override
    { return Weights_Map(); }
    bool CalculateTotalXSec(const std::string &,const bool) override
    { return false; }
    void SetLookUp(const bool) override {}
  };

  Vec4D MasslessMomentum(const double &pt,const double &y,const double &phi)
  {
    return Vec4D(pt*std::cosh(y),pt*std::cos(phi),pt*std::sin(phi),
                 pt*std::sinh(y));
  }

  struct Diphoton_Plus_Jet {
    std::unique_ptr<Flavours_Only_Process> p_proc;

    Diphoton_Plus_Jet()
    {
      InitGlobalsOnce();
      Flavour g(kf_gluon), y(kf_photon);
      p_proc.reset(new Flavours_Only_Process(Flavour_Vector{g,g,y,y,g},2));
    }

    std::unique_ptr<Selector_Base> VariableSelector(const std::string &yaml)
    {
      Selector_Key key;
      key.m_settings=Scoped_Settings("VariableSelector: {"+yaml+"}");
      key.p_proc=p_proc.get();
      return std::unique_ptr<Selector_Base>
        (Selector_Getter::GetObject("VariableSelector",key));
    }

    bool Passes(Selector_Base &selector,
                const double &pt_first_photon,const double &pt_second_photon)
    {
      Vec4D_Vector p{Vec4D(100.,0.,0.,100.),Vec4D(100.,0.,0.,-100.),
                     MasslessMomentum(pt_first_photon,0.5,0.),
                     MasslessMomentum(pt_second_photon,-0.3,2.),
                     MasslessMomentum(30.,1.,4.)};
      Selector_List sl(p_proc->Flavours(),p,p_proc->NIn());
      return selector.Trigger(sl);
    }
  };

}

TEST_CASE("Variable selector applies ranges in process order without ordering",
          "[PHASIC::Selectors::Variable_Selector]") {
  Diphoton_Plus_Jet event;
  auto selector = event.VariableSelector
    ("Variable: PT, Flavs: [22], Ranges: [[20, 1000], [18, 1000]]");
  CHECK(event.Passes(*selector,25.,19.));
  CHECK(!event.Passes(*selector,19.,25.));
  CHECK(!event.Passes(*selector,25.,17.));
  CHECK(!event.Passes(*selector,19.,19.));
}

TEST_CASE("Variable selector applies ranges to pT-ordered particles with PT_UP",
          "[PHASIC::Selectors::Variable_Selector]") {
  Diphoton_Plus_Jet event;
  auto selector = event.VariableSelector
    ("Variable: PT, Flavs: [22], Ranges: [[20, 1000], [18, 1000]],"
     " Ordering: [PT_UP]");
  CHECK(event.Passes(*selector,25.,19.));
  CHECK(event.Passes(*selector,19.,25.));
  CHECK(!event.Passes(*selector,25.,17.));
  CHECK(!event.Passes(*selector,17.,25.));
  CHECK(!event.Passes(*selector,19.,19.));
}

TEST_CASE("Variable selector with a single range and PT_UP cuts on the hardest",
          "[PHASIC::Selectors::Variable_Selector]") {
  Diphoton_Plus_Jet event;
  auto selector = event.VariableSelector
    ("Variable: PT, Flavs: [22], Ranges: [[20, 1000]], Ordering: [PT_UP]");
  CHECK(event.Passes(*selector,25.,5.));
  CHECK(event.Passes(*selector,5.,25.));
  CHECK(!event.Passes(*selector,5.,19.));
}

TEST_CASE("Variable selector rejects unknown orderings",
          "[PHASIC::Selectors::Variable_Selector]") {
  Diphoton_Plus_Jet event;
  CHECK_THROWS(event.VariableSelector
    ("Variable: PT, Flavs: [22], Ranges: [[20, 1000]],"
     " Ordering: [NOT_AN_ORDERING]"));
}
