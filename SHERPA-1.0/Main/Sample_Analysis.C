#include "Sample_Analysis.H"
#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

extern "C" {
  void pyhepc_(int&);
  void pylist_(int&);
}

Sample_Analysis::Sample_Analysis() {
  status=rpa.test.Analysis();
  cout<<" status="<<status<<endl;

}

// do some initialisation work
void Sample_Analysis::Init() {
  if (!(status)) return;

  ana = new Primitive_Analysis();
  ana->AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));

  // not (yet) included in analysis
  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me"));
  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me_nll"));
}

void Sample_Analysis::AfterME(APHYTOOLS::Blob_List * blobs) {
  // extra statistics
  obs[0]->Evaluate(*blobs,1.);
}

void Sample_Analysis::AfterPartonShower(APHYTOOLS::Blob_List * blobs) {
  msg.Events()<<" in Sample_Analysis::AfterPartonShower() "<<endl;

  int i=0;
  // fill partons FORTRAN HEPEVT
  for (Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
    if ((*bit)->Type()[0]=='F'  && (*bit)->Type()[1]=='S') {
      convert.Blob2HepEvt((*bit),1);
    }
  }
  convert.Blob2HepEvt(0,0);

  // print FORTRAN HEPEVT
  if (rpa.gen.Events()) {
    int dummy=2;           // from standard to pythia
    pyhepc_(dummy);
    dummy=1;
    pylist_(dummy);
  }

  if (!(status & 1)) return;

  // extra statistics
  obs[1]->Evaluate(*blobs,1.);

  // simple parton level analysis:
  ana->DoAnalysis(*blobs,1.);
}

void Sample_Analysis::AfterHadronization(APHYTOOLS::Blob_List * blobs) {
  if (!(status & 2)) return;


  if (rpa.gen.Events()) {
    // print FORTRAN HEPEVT
    int dummy=2;           // from standard to pythia
    pyhepc_(dummy);
    dummy=1;
    pylist_(dummy);
  }

}

void Sample_Analysis::Finish() {
  if (status== 1) {
    for (int i=0; i<obs.size();++i) {
      ana->AddObservable(obs[i]);
    }

    //    ana->FinishAnalysis("testout_sherpa_GE130c",0);
    MyStrStream s1;
    int   alf = int(1000.*rpa.consts.FixedAlphaS()+0.5);
    string salf;
    s1<<alf;
    s1>>salf;

    MyStrStream s2;
    int   yf = int(1000.*rpa.test.FactorYcut()+0.5);
    string syf;
    s2<<yf;
    s2>>syf;

    MyStrStream s3;
    int   nllf = int(100.*rpa.test.FactorNLLQ()+0.5);
    string snllf;
    s3<<nllf;
    s3>>snllf;

    string name=string("sh_GO_") + salf + string("_") + syf + string("_") + snllf;

    //    ana->FinishAnalysis("testout_sherpa_GE125g",0);
    ana->FinishAnalysis(name,0);
  }

}
