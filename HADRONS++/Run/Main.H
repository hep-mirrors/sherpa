#ifndef HADRONS_Run_Main_H
#define HADRONS_Run_Main_H

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Read_Write_Base.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Math/Random.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/ISR_Base.H"
#include "PDF/Main/Intact.H"
#include "BEAM/Main/Monochromatic.H"

#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"

#include "AHADIC++/Tools/Hadron_Init.H"

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TROOT.h"
#endif

using namespace std;
using namespace HADRONS;
using namespace ATOOLS;

void InitialiseGenerator(int argc, char *argv[]);
void InitialiseAnalysis();
Blob_List* GenerateEvent();
void AnalyseEvent(Blob_List*);
void CleanUpEvent(Blob_List*);
void FinishAnalysis();
void FinishGenerator();

static int m_analysis;

void small_sherpa_init(int argc, char *argv[])
{
  // use this for small programs, which don't use the full Sherpa framework,
  // but only Hadrons

  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Exception_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::ran = new Random(1234);
  ATOOLS::rpa = new Run_Parameter();
  ATOOLS::s_loader = new Library_Loader();

  PRINT_INFO("using small_sherpa_init");
  rpa->gen.SetVariable("PATH_PIECE","./");

  for (int i=1; i<argc;++i) {
    string par = string(argv[i]);
    string key,value;
    int equal  = par.find("=");
    if (equal!=-1) {
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
      if (key[key.length()-1]==':') {
	key.erase(key.length()-1,1);
	Read_Write_Base::AddGlobalTag(key,value);
      }
      else {
	Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
      }
    }
  }
  rpa->Init("","Run.dat",argc,argv);
  
  s_loader->LoadLibrary("SherpaSM");
  MODEL::s_model=MODEL::Model_Base::Model_Getter_Function::
    GetObject("SM",MODEL::Model_Arguments("","Model.dat",true));
  PDF::ISR_Base ** isrbases = new PDF::ISR_Base*[2];
  isrbases[0] = new PDF::Intact(Flavour(kf_e));
  isrbases[1] = new PDF::Intact(Flavour(kf_e).Bar());
  PDF::ISR_Handler_Map isrhandlers;
  isrhandlers[PDF::isr::hard_process] = new PDF::ISR_Handler(isrbases);
  isrhandlers[PDF::isr::hard_process]->SetBeam(new BEAM::Monochromatic(Flavour(kf_e),7000.,0.0,1),0);
  isrhandlers[PDF::isr::hard_process]->SetBeam(new BEAM::Monochromatic(Flavour(kf_e).Bar(),7000.,0.0,-1),1);
  double bunch_splimits[2] = {1e-10, 1.};
  isrhandlers[PDF::isr::hard_process]->Init(bunch_splimits);
  MODEL::s_model->ModelInit(isrhandlers);

  AHADIC::Hadron_Init hi;
  hi.Init();
  Default_Reader reader;
  hi.OverWriteProperties(reader);
  ATOOLS::OutputHadrons(std::cout);

  m_analysis = reader.GetValue<int>("ANALYSIS",1);

  msg->SetModifiable(true);
}


int main(int argc, char *argv[])
{
#ifdef USING__MPI
  MPI::Init(argc,argv);
#endif
  try {
#ifdef USING__ROOT
    TApplication myapp(string("HadronsApp").c_str(),&argc,&(*argv));
    gROOT->SetStyle("Plain");
#endif
    InitialiseGenerator(argc, argv);
    if(m_analysis) InitialiseAnalysis();

    for(int i=0; i<rpa->gen.NumberOfEvents(); i++) {
      Blob_List* blobs = GenerateEvent();
      if(m_analysis) AnalyseEvent(blobs);
      CleanUpEvent(blobs);

      rpa->gen.SetNumberOfGeneratedEvents(i+1);
      int steps = int(double(rpa->gen.NumberOfEvents()/10.0));
      if(steps==0) steps=1;
      if((i+1)%steps==0) msg_Info()<<i+1<<" events passed."<<std::endl;
    }
    msg_Info()<<"Generated "<<rpa->gen.NumberOfEvents()<<" events"<<endl;

    FinishGenerator();
    if(m_analysis) FinishAnalysis();
  }
  catch (Exception Sherpaexception) {
    msg_Error()<<Sherpaexception<<endl;
    terminate();
  }
  catch (std::exception stdexception) {
    cout<<"Sherpa: throws exception "
        <<stdexception.what()<<" ..."<<endl;
    terminate();
  }
#ifdef USING__MPI
  MPI::Finalize();
#endif
  return 0;
}

#ifdef USING__ROOT
TH1D* makeTH1D(string name, string title, int nbins, double xmin, double xmax,
               string xtitle, string ytitle)
{
  TH1D* hist = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->SetXTitle(xtitle.c_str());
  hist->SetYTitle(ytitle.c_str());
  hist->SetStats(kFALSE);
  return hist;
}
#endif

#endif
