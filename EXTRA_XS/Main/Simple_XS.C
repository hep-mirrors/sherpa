#ifndef EXTRA_XS_Main_Simple_XS_H
#define EXTRA_XS_Main_Simple_XS_H

#include "EXTRA_XS/Main/Process_Group.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace MODEL    { class Model_Base;   }
//namespace REMNANTS { class Remnant_Base; }

namespace EXTRAXS {

  class Simple_XS: public Process_Group, public PHASIC::ME_Generator_Base {
  private :

    void DrawLogo(std::ostream &ostr);

  public :

    // constructor
    Simple_XS();

    // destructor
    ~Simple_XS();

    // member functions
    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    Process_Base *InitializeProcess(const PHASIC::Process_Info &, bool add);
    int PerformTests();
    bool NewLibraries();

  }; // end of class Simple_XS

} // end of namespace EXTRAXS

#endif

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Model_Base.H"
//#include "REMNANTS/Main/Remnant_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "EXTRA_XS/Main/Single_Process.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

void Simple_XS::DrawLogo(std::ostream &ostr)
{
}

Simple_XS::Simple_XS(): 
  ME_Generator_Base("Internal")
{
  DrawLogo(std::cout);
}

Simple_XS::~Simple_XS() 
{
}

bool Simple_XS::Initialize(Model_Base *const model,
			   BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isrhandler)
{
  SetPSMasses();
  p_int->SetBeam(beam); 
  p_int->SetISR(isrhandler);
  return true;
}

Process_Base *Simple_XS::InitializeProcess(const Process_Info &pi, bool add)
{
  bool oneisgroup(pi.m_ii.IsGroup()||pi.m_fi.IsGroup());
  if (oneisgroup) {
    Process_Group* newxs = new Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    if (!newxs->ConstructProcesses()) {
      msg_Debugging()<<METHOD<<": group construction failed for "
		     <<newxs->Name()<<"\n";
      delete newxs;
      return NULL;
    }
    if (add) Add(newxs,1);
    newxs->SetGenerator(this);
    DEBUG_INFO("Initialized '"<<newxs->Name());
    return newxs;
  }
  else {
    Single_Process* newxs = new Single_Process();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    if (!newxs->Initialize()) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (add) Add(newxs);
    newxs->SetGenerator(this);
    DEBUG_INFO("Initialized '"<<newxs->Name());
    return newxs;
  }
}

int Simple_XS::PerformTests()
{
  return 1;
}

bool Simple_XS::NewLibraries()
{
  return false;
}

DECLARE_GETTER(Simple_XS,"Internal",
	       ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,Simple_XS>::
operator()(const ME_Generator_Key &key) const
{
  return new Simple_XS();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Simple_XS>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The internal ME generator"; 
}
