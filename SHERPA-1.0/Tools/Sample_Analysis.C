#include "Sample_Analysis.H"
#include "Shower_Observables.H"
#include "One_Particle_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"



using namespace SHERPA;
using namespace ATOOLS;
using namespace std;


int Sample_Analysis::m_njet=2;
double Sample_Analysis::m_pt_W=0.;

extern double amegic_apacic_interface_last_hard_scale;


PHard_Observable::PHard_Observable(int _type,double _xmin,double _xmax, int _nbins, 
				   std::string _name = std::string("pt_hard.dat")) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL) { name = _name; }


void PHard_Observable::Evaluate(double w)
{
  histo->Insert(Sample_Analysis::m_pt_W,w); 
}

void PHard_Observable::Evaluate(int,ATOOLS::Vec4D *,ATOOLS::Flavour *,double w)
{
  Evaluate(w);
}

void PHard_Observable::Evaluate(const ATOOLS::Particle_List &,double w)
{
  Evaluate(w);
}

void PHard_Observable::Evaluate(const ATOOLS::Blob_List &,double w)
{
  Evaluate(w);
}

Sample_Analysis::Sample_Analysis(std::string _m_path,std::string _m_file):
  m_path(_m_path),
  m_file(_m_file)
{
  status     = rpa.gen.Analysis();
  if (!(status)) { ana = 0; return; }
  ana        = new Primitive_Analysis();
  ana->AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));
  //  ana->AddObservable();
  //  ana->AddObservable(new PHard_Observable(00,0.,250.,125));
  //  ana->AddObservable(new PT_Distribution(00,0.,250.,125,1,Flavour(kf::photon)));
//   ana->AddObservable(new PT_Distribution(00,0.,250.,50,6,Flavour(kf::jet)));

  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me"));
  obs.push_back(new ME_Rate(00,1.5,7.5,6,"me_nll"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_2"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_3"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_4"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_5"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_6"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure_2"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure_3"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure_4"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure_5"));
  obs.push_back(new PHard_Observable(00,0.,200.,100,"pw_me_pure_6"));
}



void Sample_Analysis::Init() { 
  if (hepevt) nhep = 0;
}

void Sample_Analysis::AfterME(ATOOLS::Blob_List * blobs, double weight) {
  if (!(status)) return;
  obs[0]->Evaluate(*blobs,weight);

  for (Blob_Const_Iterator bit=blobs->begin();bit!=blobs->end();++bit) {
    if ((*bit)->Type().find("Signal") !=std::string::npos ) {
      m_njet= (*bit)->NOutP();
      Vec4D sum;
      for (int i=0; i<m_njet; ++i) {
	if ( (*bit)->OutParticle(i)->Flav().IsLepton()) {
	  sum+=(*bit)->OutParticle(i)->Momentum();
	}
	m_pt_W=sqrt(sqr(sum[1])+sqr(sum[2]));
      }
    }
  }
  obs[8]->Evaluate(*blobs,weight);
  obs[m_njet+7]->Evaluate(*blobs,weight);

}

void Sample_Analysis::AfterPartonShower(ATOOLS::Blob_List * blobs, double weight) {
  if (!(status & 1)) return;

  // extra statistics
  obs[1]->Evaluate(*blobs,weight);
  obs[2]->Evaluate(*blobs,weight);
  obs[m_njet+1]->Evaluate(*blobs,weight);

  
  // simple parton level analysis:
  ana->DoAnalysis(*blobs,weight);
}

void Sample_Analysis::AfterHadronization(ATOOLS::Blob_List * blobs, double weight) {
  if (!(status == 2)) return;
  ana->DoAnalysis(*blobs,weight);
}



void Sample_Analysis::Finish(std::string addpath) 
{
  if (status== 1) {
    for (size_t i=0; i<obs.size();++i) {
      ana->AddObservable(obs[i]);
    }

    MyStrStream s1;
    int   alf = int(1000.*rpa.gen.ScalarFunction(string("alpha_S"))+0.5);
    string salf;
    s1<<alf;
    s1>>salf;

    ATOOLS::Data_Read *dataread = new Data_Read(m_path+m_file);
    std::string outputpath=dataread->GetValue<std::string>("ANALYSIS_OUTPUT",std::string("output"));
    delete dataread;
    msg.Out()<<" FinishAnalysis("<<outputpath+addpath<<");"<<endl;
    ana->FinishAnalysis(outputpath+addpath,0);
  }

}


Sample_Analysis::~Sample_Analysis() {
  if (ana)                   { delete ana;     ana     = 0; }
}
