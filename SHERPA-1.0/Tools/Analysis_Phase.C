#include "Analysis_Phase.H"
#include "Shower_Observables.H"
#include "One_Particle_Observables.H"
#include "Two_Particle_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"



using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Analysis_Phase::Analysis_Phase() :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(std::string("")), p_analysis(NULL) { }

Analysis_Phase::Analysis_Phase(std::string _btype) :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(_btype), p_analysis(NULL) 
{
  m_type        = string("Perturbative");
  m_status      = 0;
  p_analysis    = new Primitive_Analysis(m_btype);
  Flavour flav1 = Flavour(kf::nue), flav2 = flav1.Bar();
  p_analysis->AddObservable(new One_Particle_PT(flav1,00,0.,200.,100));
  p_analysis->AddObservable(new One_Particle_Eta(flav1,00,-3.,3.,30));
  p_analysis->AddObservable(new Two_Particle_Mass(flav1,flav2,00,0.,200.,100));
  p_analysis->AddObservable(new Two_Particle_PT(flav1,flav2,00,0.,200.,100));
  p_analysis->SetBlobType(_btype);
}

Analysis_Phase::Analysis_Phase(Primitive_Analysis * _ana,std::string _btype) :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(_btype), p_analysis(_ana) 
{
  m_type      = string("Perturbative");
  m_status    = 0;
}





bool  Analysis_Phase::Treat(ATOOLS::Blob_List * _blist, double & _weight) 
{
  p_analysis->DoAnalysis(*_blist,1.);
  return 0;
}


void  Analysis_Phase::CleanUp() 
{
  m_status=0;
}

void Analysis_Phase::Finish(std::string _dirname)
{
  std::cout<<"In Analysis_Phase::Finish : "<<_dirname<<std::endl;
  p_analysis->FinishAnalysis(_dirname);
}
