#include "Sample_Analysis.H"
#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"



namespace SHERPA {

  class PHard_Observable : public ATOOLS::Primitive_Observable_Base {  
  public:
    PHard_Observable(int _type,double _xmin,double _xmax,int _nbins, std::string _name);
    void Evaluate(double); 
    void Evaluate(int,ATOOLS::Vec4D *,ATOOLS::Flavour *,double w);
    void Evaluate(const ATOOLS::Particle_List &,double w);
    void Evaluate(const ATOOLS::Blob_List &,double w); 
  };

  extern double amegic_apacic_interface_last_hard_scale;
}

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

int Sample_Analysis::m_njet=2;
double Sample_Analysis::m_pt_W=0.;

PHard_Observable::PHard_Observable(int _type,double _xmin,double _xmax,int _nbins, std::string _name)
{
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins;
  name  = std::string("pt_hard.dat");
  name  = _name+std::string(".dat");
  histo = new ATOOLS::Histogram(type,xmin,xmax,nbins);
}

void PHard_Observable::Evaluate(double w)
{
  //  histo->Insert(sqrt(amegic_apacic_interface_last_hard_scale),w); 
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


void Sample_Analysis::Init() { 
  if (hepevt) nhep = 0;
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
    m_type      = std::string("Hadronization");
    msg.Tracking()<<" Init "<<m_name<<endl;
    break;
  default:
    msg.Out()<<" Unknow Analysis_Phase "<<m_mode<<endl;
  }
  if (!ana) m_mode=0;

}





bool  Analysis_Phase::Treat(ATOOLS::Blob_List * bl, double & weight) 
{
  switch (m_mode) {
  case 1:
    // should be if "signal process" has status 1! beware of MI
    if (bl->size()==1 && bl->back()->Status()==1) {

      p_analysis->AfterME(bl,weight);
      m_status = 1;
    }
    break;
  case 2:
    if (m_status) return 0;
    // should look for a shower blobs! 
    if (bl->size()>2) {
      p_analysis->AfterPartonShower(bl,weight);
      m_status = 1;
    }
    break;
  case 3:
    if (m_status) return 0;
    if (bl->size()>2) {
      p_analysis->AfterHadronization(bl,weight);
      m_status = 1;
    }
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

