#include "Analysis_Phase.H"
#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"
#include <vector>

#include "One_Particle_Observables.H"
#include "Two_Particle_Observables.H"
#include "Four_Particle_Observables.H"

#include "Shower_Observables.H"
#include "Primitive_Detector.H"
#include "Jet_Observables.H"
#include "Final_Selector.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif
using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Analysis_Phase::Analysis_Phase() :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(std::string("")), p_analysis(NULL) { }

Analysis_Phase::Analysis_Phase(const std::string & btype,const std::string & path,const std::string & file) :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(btype), p_analysis(NULL), m_path(path), m_file(file) 
{
  //  m_type        = string("Perturbative");
  m_type        = string("Hadronization");
  m_status      = 0;
  p_analysis    = new Primitive_Analysis(m_btype,ANALYSIS::fill_all|ANALYSIS::splitt_jetseeds|
					 ANALYSIS::splitt_phase|ANALYSIS::do_me|ANALYSIS::do_shower);
  //ANALYSIS::splitt_phase|ANALYSIS::do_me);
  // --- Andreas Analysis ---

  Final_Selector * fsel =  new Final_Selector("FinalState","KtJetsNLeptons");
  // pure teilchen eigenschaften
  Final_Selector_Data fd;
  //LHC analysis
  
  //fd.rmin   =0.4;
  fd.rmin   =1.;
  fd.pt_min =20.;
  fd.eta_min=-4.5;
  fd.eta_max=+4.5;
  
  //tevatron analysis
  /*
  fd.rmin   =1.;
  fd.pt_min =15.;
  fd.eta_min=-2.;
  fd.eta_max=+2.;
  */
  fsel->AddSelector(Flavour(kf::jet),fd);
  
  //LHC analysis
  Final_Selector_Data fcorrlj;
  fcorrlj.rmin = 0.4;
  fsel->AddSelector(Flavour(kf::jet),Flavour(kf::lepton),fcorrlj);
  
  Final_Selector_Data fcorrll;
  fcorrll.rmin = 0.2;
  fsel->AddSelector(Flavour(kf::lepton),Flavour(kf::lepton),fcorrll);
  
  p_analysis->AddObservable(fsel);

  // all jets
  // bedingung 2 lepton
  fsel = new Final_Selector("KtJetsNLeptons","KtJets");
  fsel->AddKeepFlavour(Flavour(kf::jet));
  p_analysis->AddObservable(fsel);

  // all lepton
  fsel = new Final_Selector("KtJetsNLeptons","Leptons");
  fsel->AddKeepFlavour(Flavour(kf::lepton));
  p_analysis->AddObservable(fsel);
 

  //p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,0,"FinalState"));
  //p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,0,"KtJets"));
  /*
  p_analysis->AddObservable(new Two_Particle_PT(Flavour(kf::mu),Flavour(kf::mu).Bar(),
						00,0.,150.,150,"Zpt"));
  */
  /*
  p_analysis->AddObservable(new Two_Particle_Eta(Flavour(kf::mu),Flavour(kf::mu).Bar(),
						00,-5.,5.,50,"Zeta"));
  */
  //exclusive analysis
  //  p_analysis->AddObservable(new Jet_PT_Distribution(2,15.,200.,37,1,5,"KtJets"));
  //inclusive one
  // p_analysis->AddObservable(new Jet_PT_Distribution(1,15.,200.,37,1,5,"KtJets"));
  
  //LHC
  p_analysis->AddObservable(new Jet_PT_Distribution(1,15.,205.,38,1,5,"KtJets"));
  //tevatron 2jets
  //p_analysis->AddObservable(new Jet_PT_Distribution(2,15.,205.,38,2,5,"KtJets"));
  
  //p_analysis->AddObservable(new Jet_PT_Distribution(1,17.5,142.5,25,0,5,"KtJets"));
  
  //p_analysis->AddObservable(new Jet_Eta_Distribution(1,-4.75,4.75,19,1,5,"KtJets"));
}

Analysis_Phase::Analysis_Phase(Primitive_Analysis * ana,const std::string  & btype) :
  Event_Phase_Handler(std::string("Analysis")), 
  m_btype(btype), p_analysis(ana) 
{
  m_type      = string("Perturbative");
  m_status    = 0;
}





bool  Analysis_Phase::Treat(ATOOLS::Blob_List * blist, double & weight) 
{
//   cout<<" Analysis:: "<<endl;
//   cout<<*blist<<endl;
    

  PROFILE_LOCAL("Analysis_Phase::Treat");
  p_analysis->DoAnalysis(blist,weight);
  return 0;
}


void  Analysis_Phase::CleanUp() 
{
  m_status=0;
}

void Analysis_Phase::Finish(const std::string & dirname)
{
  ATOOLS::Data_Read *dataread = new Data_Read(m_path+m_file);
  std::string outputpath=dataread->GetValue<std::string>("ANALYSIS_OUTPUT",std::string("output"));
  delete dataread;
  msg.Out()<<" FinishAnalysis("<<outputpath<<");"<<endl;
  p_analysis->FinishAnalysis(outputpath);
}
