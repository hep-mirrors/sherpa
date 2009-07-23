#ifndef HADRONS_Run_Main_H
#define HADRONS_Run_Main_H

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Read_Write_Base.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "MODEL/Main/Standard_Model.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"

#include "HADRONS++/Main/Hadrons.H"
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


void small_sherpa_init(int argc, char *argv[])
{
  // use this for small programs, which don't use the full Sherpa framework,
  // but only Hadrons

  PRINT_INFO("using small_sherpa_init");
  std::vector<std::string> searchpaths;
  searchpaths.push_back(rpa.gen.Variable("SHERPA_RUN_PATH")+"/");
  My_Out_File::SetSearchPaths(searchpaths);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH")+"/");
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH"));
  searchpaths.push_back(SHERPA_SHARE_PATH+std::string("/"));
  searchpaths.push_back(SHERPA_SHARE_PATH);
  My_In_File::SetSearchPaths(searchpaths);
  rpa.gen.SetVariable("PATH_PIECE","./");

  std::vector<std::string> names(4);
  names[0]="Decaydata";
  names[1]="Run.dat";
  My_In_File::SetNoComplains(names);

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
  ATOOLS::s_loader = new Library_Loader();
  rpa.Init("","Run.dat",argc,argv);
  MODEL::s_model = new MODEL::Standard_Model("","Model.dat",true);
  MODEL::s_model->InitializeInteractionModel();
  
  AHADIC::Hadron_Init hi;
  hi.Init();
  Data_Reader overwrite;
  hi.OverWriteProperties(overwrite);

  msg->SetModifiable(true);
}


int main(int argc, char *argv[])
{
  PRINT_INFO("starting main");
  ATOOLS::exh->Init();
  set_terminate(ATOOLS::Terminate);
  set_unexpected(ATOOLS::Terminate);
  signal(SIGSEGV,ATOOLS::SignalHandler);
  signal(SIGINT,ATOOLS::SignalHandler);
  signal(SIGBUS,ATOOLS::SignalHandler);
  signal(SIGFPE,ATOOLS::SignalHandler);
  signal(SIGABRT,ATOOLS::SignalHandler);
  signal(SIGTERM,ATOOLS::SignalHandler);
  signal(SIGXCPU,ATOOLS::SignalHandler);
  try {
#ifdef USING__ROOT
    TApplication myapp(string("HadronsApp").c_str(),&argc,&(*argv));
    gROOT->SetStyle("Plain");
#endif
    PRINT_INFO("starting generator");
    InitialiseGenerator(argc, argv);
    PRINT_INFO("starting analysis");
    if(rpa.gen.Analysis()) InitialiseAnalysis();

    for(int i=0; i<rpa.gen.NumberOfEvents(); i++) {
      Blob_List* blobs = GenerateEvent();
      if(rpa.gen.Analysis()) AnalyseEvent(blobs);
      CleanUpEvent(blobs);

      rpa.gen.SetNumberOfDicedEvents(i+1);
      int steps = int(double(rpa.gen.NumberOfEvents()/10.0));
      if(steps==0) steps=1;
      if((i+1)%steps==0) msg_Info()<<i+1<<" events passed."<<std::endl;
    }
    msg_Info()<<"Generated "<<rpa.gen.NumberOfEvents()<<" events"<<endl;

    FinishGenerator();
    if(rpa.gen.Analysis()) FinishAnalysis();
  }
  catch (Exception Sherpaexception) {
    Sherpaexception.UpdateLogFile();
    msg_Error()<<Sherpaexception<<endl;
    terminate();
  }
  catch (std::exception stdexception) {
    cout<<"Sherpa: throws exception "
        <<stdexception.what()<<" ..."<<endl;
    terminate();
  }
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
