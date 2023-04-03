#ifndef NEUTRINOS_Main_Tests_H
#define NEUTRINOS_Main_Tests_H

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Math/Random.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/ISR_Base.H"
#include "PDF/Main/Intact.H"
#include "BEAM/Spectra/Monochromatic.H"
#include "SHERPA/SoftPhysics/Hadron_Init.H"
#include "SHERPA/Initialization/Initialization_Handler.H"


#include "NEUTRINOS++/Current_Library/Current_ME.H"
#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"

void small_sherpa_init(int argc, char *argv[])
{
  // use this for small programs, which don't use the full Sherpa framework,
  // but only Hadrons
  ATOOLS::mpi      = new ATOOLS::My_MPI();
  ATOOLS::exh      = new ATOOLS::Terminator_Object_Handler();
  ATOOLS::msg      = new ATOOLS::Message();
  ATOOLS::Settings::InitializeMainSettings(argc, argv);
  ATOOLS::ran      = new ATOOLS::Random(1234);
  ATOOLS::rpa      = new ATOOLS::Run_Parameter();
  ATOOLS::s_loader = new ATOOLS::Library_Loader();

  PRINT_INFO("using small_sherpa_init");

  ATOOLS::rpa->Init();

  ATOOLS::s_loader->LoadLibrary("SherpaSM");
  SHERPA::Initialization_Handler init;
  MODEL::s_model=MODEL::Model_Base::Model_Getter_Function::
    GetObject("SM",MODEL::Model_Arguments(true));
  PDF::ISR_Base ** isrbases = new PDF::ISR_Base*[2];
  isrbases[0] = new PDF::Intact(ATOOLS::Flavour(kf_e));
  isrbases[1] = new PDF::Intact(ATOOLS::Flavour(kf_e).Bar());
  PDF::ISR_Handler_Map isrhandlers;
  isrhandlers[PDF::isr::hard_process] = new PDF::ISR_Handler(isrbases);
  isrhandlers[PDF::isr::hard_process]->SetBeam(new BEAM::Monochromatic(ATOOLS::Flavour(kf_e),7000.,0.0,1),0);
  isrhandlers[PDF::isr::hard_process]->SetBeam(new BEAM::Monochromatic(ATOOLS::Flavour(kf_e).Bar(),7000.,0.0,-1),1);
  MODEL::s_model->ModelInit(isrhandlers);

  SHERPA::Hadron_Init().Init();
  //ATOOLS::OutputHadrons(std::cout);
  ATOOLS::msg->SetModifiable(true);
}


namespace NEUTRINOS {
  class Tests {
  private:
    ATOOLS::Flavour_Vector     m_flavs;   
    std::vector<int>           m_indices; 
    std::vector<ATOOLS::Vec4D> m_momenta; 

    Current_ME           * p_ME;
    Scatter_Current_Base * p_JLL, * p_JNN;
    
    void InitVectors(const size_t n=4);
    void MakeSimpleScatterKinematics(const double & E,const double & costheta,const double & phi);
    void ep_ep();
  public:
    Tests(int argc,char* argv[]);
    ~Tests() {};
  };
}
#endif
