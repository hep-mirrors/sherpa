#include "Sample_Analysis.H"
#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

/*
  extern "C" {
  void pyhepc_(int&);
  void pylist_(int&);
  }
*/


void Sample_Analysis::Init() { 
  if (hepevt) nhep = 0;
}

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
  if (!(status)) { ana = 0; return; }
  ana        = new Primitive_Analysis();
  ana->AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));
//   ana->AddObservable(new PT_Distribution(00,0.,250.,50,1,Flavour(kf::W)));
//   ana->AddObservable(new PT_Distribution(00,0.,250.,50,6,Flavour(kf::jet)));

  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me"));
  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me_nll"));
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
  // extra statistics
  //  cout<<" in Sample_Analysis::AfterME with "<<blobs->size()<<" ("<<status<<")"<<endl;
  if (!(status)) return;
  obs[0]->Evaluate(*blobs,1.);
}

void Sample_Analysis::AfterPartonShower(APHYTOOLS::Blob_List * blobs) {
  // fill partons into FORTRAN HEPEVT and print it 
  //  cout<<" in Sample_Analysis::AfterPartonShower with "<<blobs->size()<<" ("<<status<<")"<<endl;
  if (hepevt) {
    convert->Blobs2HepEvt(blobs,std::string("ISR"),nhep);
    convert->Blobs2HepEvt(blobs,std::string("FSR"),nhep);

    msg.Events()<<" HepEVT after Parton Shower"<<endl;
    //pyhepc_(hepevtmode);    // from standard to pythia
    //pylist_(hepevtmode);
  }

  if (!(status & 1)) return;

  // extra statistics
  obs[1]->Evaluate(*blobs,1.);

  // simple parton level analysis:
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

  if (!(status & 2)) return;
  //  ana->DoAnalysis(*blobs,1.);
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
    int   yf = int(100.*rpa.test.FactorYcut()+0.5);
    string syf;
    s2<<yf;
    s2>>syf;

    MyStrStream s3;
    int   nllf = int(100.*rpa.test.FactorNLLQ()+0.5);
    string snllf;
    s3<<nllf;
    s3>>snllf;

    string name=string("sh_TR_") + salf + string("_") + syf + string("_") + snllf;

    msg.Out()<<" FinishAnalysis("<<name<<");"<<endl;
    //    ana->FinishAnalysis("testout_sherpa_GE125g",0);
    ana->FinishAnalysis(name,0);
  }

}


Sample_Analysis::~Sample_Analysis() {
  if (ana)                   { delete ana;     ana     = 0; }
  if ((convert) && (InitIO)) { delete convert; convert = 0; } 
}


// ============================================================


Analysis_Phase::Analysis_Phase(Sample_Analysis * ana, int mode) :
  p_analysis(ana), m_mode(mode) 
{
  m_type      = string("Perturbative");
  m_status    = 0;

  if (!ana) m_mode=0;

  switch (m_mode) {
  case 1:
    m_name      = string("Analysis_Phase Matrix Elements");
    msg.Tracking()<<" Init "<<m_name<<endl;
    break;
  case 2:
    m_name      = string("Analysis_Phase Parton Shower");
    msg.Tracking()<<" Init "<<m_name<<endl;
    break;
  case 3:
    m_name      = string("Analysis_Phase Hadronization");
    msg.Tracking()<<" Init "<<m_name<<endl;
    break;
  default:
    msg.Out()<<" Unknow Analysis_Phase "<<m_mode<<endl;
  }
  if (!ana) m_mode=0;

}





bool  Analysis_Phase::Treat(APHYTOOLS::Blob_List * bl) 
{
  switch (m_mode) {
  case 1:
    if (bl->size()==1 && bl->back()->Status()==1) {
      p_analysis->AfterME(bl);
      m_status = 1;
    }
    break;
  case 2:
    if (m_status) return 0;
    if (bl->size()>1) {
      p_analysis->AfterPartonShower(bl);
      m_status = 1;
    }
    break;
  case 3:
    if (m_status) return 0;
    p_analysis->AfterHadronization(bl);
    m_status = 1;
    break;
  default:
    msg.Out()<<" Unknow Analysis_Phase "<<m_mode<<endl;
  }

  return 0; // analysis does not create any new blobs.
}

void  Analysis_Phase::CleanUp() 
{
  m_status=0;
}

