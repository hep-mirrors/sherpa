#include "Sample_Analysis.H"
#include "Shower_Observables.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

/*
  extern "C" {
  void pyhepc_(int&);
  void pylist_(int&);
  }
*/

Sample_Analysis::Sample_Analysis(IO_HepEvt *& _convert,bool _hepevt) : 
  convert(_convert), hepevt(_hepevt)
{
  InitIO       = 0;
  if (hepevt) {
    hepevtmode = 2;
    if (convert==0) {
      convert  = new IO_HepEvt();
      _convert = convert;
      InitIO   = 1;
    }
  } 

  status     = rpa.test.Analysis();
  status     = 0;
  if (!(status)) { ana = 0; return; }
  ana        = new Primitive_Analysis();
  ana->AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));
}

Sample_Analysis::~Sample_Analysis() {
  if (ana)                   { delete ana;     ana     = 0; }
  if ((convert) && (InitIO)) { delete convert; convert = 0; } 
}

void Sample_Analysis::Init() { 
  if (hepevt) nhep = 0;
}


void Sample_Analysis::AfterME(APHYTOOLS::Blob_List * blobs) {
  // fill partons into FORTRAN HEPEVT and print it 
  if (hepevt) {
    convert->Blobs2HepEvt(blobs,std::string("Beam"),nhep);
    convert->Blobs2HepEvt(blobs,std::string("ISR"),nhep);
    convert->Blobs2HepEvt(blobs,std::string("Hard ME"),nhep);

    msg.Events()<<" HepEVT after Matrix Element"<<endl;
    //pyhepc_(hepevtmode);    // from standard to pythia
    //pylist_(hepevtmode);
  }

  if (!(status & 1)) return;
  ana->DoAnalysis(*blobs,1.);
}


void Sample_Analysis::AfterPartonShower(APHYTOOLS::Blob_List * blobs) {
  // fill partons into FORTRAN HEPEVT and print it 
  if (hepevt) {
    convert->Blobs2HepEvt(blobs,std::string("ISR"),nhep);
    convert->Blobs2HepEvt(blobs,std::string("FSR"),nhep);

    msg.Events()<<" HepEVT after Parton Shower"<<endl;
    //pyhepc_(hepevtmode);    // from standard to pythia
    //pylist_(hepevtmode);
  }

  if (!(status & 1)) return;
  ana->DoAnalysis(*blobs,1.);
}


void Sample_Analysis::AfterHadronization(APHYTOOLS::Blob_List * blobs) {
  // fill partons into FORTRAN HEPEVT and print it 
  if (hepevt) {
    convert->Blobs2HepEvt(blobs,std::string("ISR"),nhep);
    convert->Blobs2HepEvt(blobs,std::string("FSR"),nhep);

    msg.Events()<<" HepEVT after Parton Shower"<<endl;
    //pyhepc_(hepevtmode);    // from standard to pythia
    //pylist_(hepevtmode);
  }

  if (!(status & 1)) return;
  ana->DoAnalysis(*blobs,1.);
}



void Sample_Analysis::Finish() {
  if (status== 1) ana->FinishAnalysis("parton_testout",0);
}
