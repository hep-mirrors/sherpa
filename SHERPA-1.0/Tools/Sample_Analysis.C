#include "Sample_Analysis.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"


#include "Shower_Observables.H"
#include "Primitive_Detector.H"
#include "One_Particle_Observables.H"


namespace SHERPA {

  class PHard_Observable : public ATOOLS::Primitive_Observable_Base {  
  public:
    PHard_Observable(int type,double xmin,double xmax,int nbins, const std::string & name);
    void Evaluate(double,int); 
    void Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double w, int);
    void Evaluate(const ATOOLS::Particle_List &,double w, int);
    void Evaluate(const ATOOLS::Blob_List &,double w, int); 
    Primitive_Observable_Base * Copy() const;
  };

  extern double amegic_apacic_interface_last_hard_scale;
}

using namespace SHERPA;
using namespace ATOOLS;
//using namespace ATOOLS::ANALYSIS;
using namespace std;


int Sample_Analysis::m_njet=2;
double Sample_Analysis::m_pt_W=0.;

PHard_Observable::PHard_Observable(int type,double xmin,double xmax,int nbins, const std::string & name)
{
  m_splitt_flag=false;
  m_type = type; m_xmin = xmin; m_xmax = xmax; m_nbins = nbins;
  m_name  = std::string("pt_hard.dat");
  m_name  = name+std::string(".dat");
  p_histo = new ATOOLS::Histogram(m_type,m_xmin,m_xmax,m_nbins);
}



void PHard_Observable::Evaluate(double w, int ncount)
{
  //  histo->Insert(sqrt(amegic_apacic_interface_last_hard_scale),w); 
  p_histo->Insert(Sample_Analysis::m_pt_W,w,ncount); 
}

void PHard_Observable::Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double w, int ncount)
{
  Evaluate(w,ncount);
}

void PHard_Observable::Evaluate(const ATOOLS::Particle_List &,double w, int ncount)
{
  Evaluate(w,ncount);
}

void PHard_Observable::Evaluate(const ATOOLS::Blob_List &,double w, int ncount)
{
  Evaluate(w,ncount);
}

Primitive_Observable_Base * PHard_Observable::Copy() const {
  std::cerr<<" ERROR in PHard_Observable::Copy() : not splittable "<<std::endl;
  return 0;
}


void Sample_Analysis::Init() { 
  if (hepevt) nhep = 0;
}

Sample_Analysis::Sample_Analysis(const std::string & path,const std::string & file):
  m_path(path),
  m_file(file)
{
  status     = rpa.gen.Analysis();
  if (!(status)) { ana = 0; return; }
  ana        = new Primitive_Analysis("Sample_Analysis",ANALYSIS::splitt_all|ANALYSIS::fill_all);
  //  ana->AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));
  ana->AddObservable(new Primitive_Detector(0.5,8.,0,"ConeJets(.5)"));
  ana->AddObservable(new Jetrates(11,1.e-6,1.,80,0));
  ana->AddObservable(new Jetrates(11,1.e-6,1.,80,0,"ConeJets(.5)"));
  ana->AddObservable(new Multiplicity(00,-0.5,50.5,51,0,"FinalState"));
  ana->AddObservable(new Multiplicity(00,-0.5,50.5,51,0,"ConeJets(.5)"));
  ana->AddObservable(new Multiplicity(00,-0.5,50.5,51,0,"KtJets"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,1,Flavour(kf::W)));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,6,Flavour(kf::jet)));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,0,4,Flavour(kf::jet),"ConeJets(.5)"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,1,4,Flavour(kf::jet),"ConeJets(.5)"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,2,4,Flavour(kf::jet),"ConeJets(.5)"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,3,4,Flavour(kf::jet),"ConeJets(.5)"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,4,4,Flavour(kf::jet),"ConeJets(.5)"));
  ana->AddObservable(new PT_Distribution(00,0.,200.,100,6,Flavour(kf::jet),"KtJets"));


  //  ana->AddObservable(new PT_Distribution(00,0.,250.,50,6,Flavour(kf::Z)));

  //  ana->AddObservable(new PHard_Observable(00,0.,250.,125));
  //  ana->AddObservable(new PT_Distribution(00,0.,250.,125,1,Flavour(kf::photon)));

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
  ana->DoAnalysis(blobs,weight);
}

void Sample_Analysis::AfterHadronization(ATOOLS::Blob_List * blobs, double weight) {
  if (!(status == 2)) return;
  ana->DoAnalysis(blobs,weight);
}



void Sample_Analysis::Finish(std::string addpath) 
{
  if (status>= 1) {
    Primitive_Analysis * sana=new Primitive_Analysis("ME_Analysis",ANALYSIS::fill_all);
    for (size_t i=0; i<obs.size();++i) {
      sana->AddObservable(obs[i]);
    }
    ana->AddSubAnalysis("ME",sana);

    MyStrStream s1;
    int   alf = int(1000.*rpa.gen.ScalarFunction(string("alpha_S"))+0.5);
    string salf;
    s1<<alf;
    s1>>salf;

    ATOOLS::Data_Read *dataread = new Data_Read(m_path+m_file);
    std::string outputpath=dataread->GetValue<std::string>("ANALYSIS_OUTPUT",std::string("output"));
    delete dataread;
    msg.Out()<<" FinishAnalysis("<<outputpath+addpath<<");"<<endl;
    ana->FinishAnalysis(outputpath+addpath);
  }

}


Sample_Analysis::~Sample_Analysis() {
  if (ana)   { delete ana;     ana     = 0; }
}
