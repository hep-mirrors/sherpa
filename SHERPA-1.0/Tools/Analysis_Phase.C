#include "Analysis_Phase.H"
#include "Shower_Observables.H"
#include "One_Particle_Observables.H"
#include "Two_Particle_Observables.H"
#include "Four_Particle_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"
#include <vector>


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
  Flavour flav1 = Flavour(kf::u), flav2 = Flavour(kf::d).Bar(),
    flav3 = Flavour(kf::u).Bar(), flav4 = Flavour(kf::d);
  std::vector<Flavour> flavs;
  flavs.push_back(flav1);
  flavs.push_back(flav2);
  flavs.push_back(flav3);
  flavs.push_back(flav4);
  p_analysis->AddObservable(new Two_Particle_Mass(flav1,flav2,00,0.,200.,100,"Mass"));
  p_analysis->AddObservable(new Two_Particle_PT(flav1,flav2,00,0.,200.,100,"PT"));
  p_analysis->AddObservable(new Two_Particle_Mass(flav3,flav4,00,0.,200.,100,"Mass"));
  p_analysis->AddObservable(new Four_Particle_PlaneAngle(flavs,00,-1.,1.,100,"PlaneAngle"));
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
  p_analysis->FinishAnalysis(_dirname);
}
