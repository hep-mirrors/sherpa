#include "Apacic.H"

#include "ISR_Handler.H"
#include "Beam_Handler.H"

#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Tree.H"
#include "Knot.H"

#include "Run_Parameter.H"
#include "Parton_List.H"
#include "Blob_List.H"
#include "Flavour.H"
#include "Random.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"


using namespace APACIC;
using namespace ISR;
using namespace BEAM;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;


Apacic::Apacic(bool _t) : test(_t) { isr = 0; beam = 0; };

void Apacic::Init() {
  ParticleInit("./");
  rpa.Init("./");

  if (!as)   as     = new Running_AlphaS();
  if (!aqed) aqed   = new Running_AlphaQED();

  seldata = new Selector_Data(rpa.GetPath());

  Data_Read dr(rpa.GetPath()+string("/ISR.dat"));

  int * isrtypes  = new int[2];
  isrtypes[0]     = dr.GetValue<ISR_Type::code>("ISR1");
  isrtypes[1]     = dr.GetValue<ISR_Type::code>("ISR2");

  double s          = sqr(AORGTOOLS::rpa.gen.Ecms());
  double * splimits = new double[2];
  splimits[0]       = s*sqr(dr.GetValue<double>("SMIN"));
  splimits[1]       = s*sqr(dr.GetValue<double>("SMAX"));

  Flavour * bunches, * beams, * parton;

  bunches      = new Flavour[2];
  beams        = new Flavour[2];
  parton       = new Flavour[2];
  for (int i=0;i<2;i++) {
    bunches[i] = Flavour(kf::none);
    beams[i]   = Flavour(kf::none);
    parton[i]  = Flavour(kf::none);
  }

  bunches[0] = beams[0] = rpa.gen.Beam1();// bunches[] is just for testing ..
  bunches[1] = beams[1] = rpa.gen.Beam2();

  if (beams[0]==Flavour(kf::e) || beams[0]==Flavour(kf::e).Bar()) {
    parton[0] = beams[0];
    parton[1] = beams[1];
  }
  else if (beams[0]==Flavour(kf::p_plus) || beams[0]==Flavour(kf::p_plus).Bar()) {
    parton[0] = beams[0];
    parton[1] = beams[1];
    cout<<" This is just a Testrun ! "<<endl; 
  }
  else {
    cout<<" partons have to be set in Apacic.C "<<endl;
    cout<<" try e+ e- and start again "<<endl;
    exit(1);
  }

  int * beamtypes  = new int[2];
  beamtypes[0] = 0;
  beamtypes[1] = 0;
  double * pol = new double[2]; 
  pol[0] = 0;
  pol[1] = 0;

  if (!beam) beam = new Beam_Handler(beamtypes,bunches,pol,beams,pol,splimits);
  if (!isr)  isr  = new ISR_Handler(isrtypes,beams,parton,splimits);


  ini = 1;
  fin = 1;
  had = 0;
  
  if (test) {
    partons            = 0;
    blobs              = 0;
    me_ps_interface    = 0;
    AP_proc            = 0;

    TestShower();
  }
  else {
    partons            = new Parton_List;
    blobs              = new Blob_List;
    me_ps_interface    = new Hard_Interface(isr,2,0);
  }
}

Apacic::~Apacic() {
  if (!test) {
    CleanUpEvent();
    if (partons) {
      if (!(partons->empty())) partons->erase(partons->begin(),partons->end()); 
      delete partons;
    }
    if (blobs) {
      if (!(blobs->empty())) blobs->erase(blobs->begin(),blobs->end()); 
      delete blobs;
    }
    if (me_ps_interface) { delete me_ps_interface; me_ps_interface = 0; }
    if (AP_proc)            delete AP_proc;
  }
  if (isr)      delete isr;
  if (beam)     delete beam;
  msg.Tracking()<<"APACIC is deleting Running_AlphaS ..."<<endl;
  if (as)       delete as;
  if (aqed)     delete aqed;
  msg.Tracking()<<"Out Apacic::~Apacic()."<<std::endl;
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
}


void Apacic::TestShower()
{  
  bool decision;
  msg.Out()<<"Test Initial or Final State Shower ? (Final = 1, Initial = 0)"<<std::endl;
  std::cin>>decision;
  
  Tree *  tree_fin = new Tree();     
  Final_State_Shower * FinShower  = new Final_State_Shower();

  if (decision) {
    rpa.gen.SetBeam1(Flavour(kf::code(11)));
    rpa.gen.SetBeam2(Flavour(kf::code(11)).Bar());
    msg.Tracking()<<"Set beams on e- e+."<<std::endl;

    FinShower->TestShower(tree_fin);
    msg.Tracking()<<"Test of Final_State_Shower complete."<<std::endl;
  }
  else {
    //    msg.Tracking()<<"ISR_Handler : "<<isr->Type()<<std::endl;
    Tree ** tree_ini = new Tree*[2];     
    for (int i=0;i<2;i++) tree_ini[i] = new Tree();
    Initial_State_Shower * InShower = new Initial_State_Shower(isr,FinShower);

    InShower->TestShower(tree_ini);
    msg.Tracking()<<"Test of Initial_State_Shower complete."<<std::endl;
    if (tree_ini) {
      for (int i=0;i<2;i++) delete tree_ini[i];
      delete [] tree_ini;
    }
    if (InShower)  delete InShower;
  }
  
  if (tree_fin)  delete tree_fin;
  if (FinShower) delete FinShower;
}


/*--------------------------------------------------------

  Stuff to mimick MOCAIC features.

  --------------------------------------------------------*/

void Apacic::CrossSections() {
  bool success;
  AP_proc = new Hard_Processes(seldata,isr,beam,success);  
  if (success) {
    if (AP_proc->PrepareCalculation()) {
      AP_proc->CalculateCrossSections();
    }
  }
}

void Apacic::GenerateEvents() {
  msg.Tracking()<<"Start generating events : "<<std::endl
		<<"Use "<<AP_proc->Name()
		<<" for the hard processes in event generation."<<std::endl;
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) {
    msg.Events()<<"##########################################################"<<std::endl
		<<"#############################"<<n<<"th event #############"<<std::endl
		<<"##########################################################"<<std::endl;
    if (AP_proc->OneEvent()) {
      CleanUpEvent();
      
      APHYTOOLS::Blob * blob;
      blob = new APHYTOOLS::Blob();
      FillBlob(blob);
      blob->SetId(blobs->size());
      blobs->push_back(blob);

      if (!(me_ps_interface->ReduceBlob(blob))) {
	msg.Error()<<"Reduction of hard blob did not work out !"<<std::endl;
	continue;
      }
      if (!(me_ps_interface->PerformShowers(ini,fin))) {
	msg.Error()<<"Performing the showers did not work out !"<<std::endl;
	continue;
      } 
      me_ps_interface->ExtractPartons(ini,fin,blobs,partons);
      /*
	if (had) {
	part_had_interface->HadronsToPartons(blobs,partons);
	part_had_interface->ExtractSinglets(partons);
	}
      */
      //      msg.Events()<<(*blobs)<<std::endl;
    }
  }
  me_ps_interface->FinalStats();
}

void Apacic::CleanUpEvent() {
  msg.Debugging()<<"In Apacic::CleanUpEvent() :"<<std::endl;
  if (!blobs->empty()) {
    for (APHYTOOLS::Blob_Iterator bit=blobs->begin(); 
	 bit!=blobs->end(); ++bit) {
      if (rpa.gen.Debugging()) msg.Tracking()<<" Delete :"<<(*bit)<<std::endl;
      (*bit)->DeleteOwnedPartons();
    }
    blobs->erase(blobs->begin(),blobs->end());
    msg.Debugging()<<"   Erased and deleted blobs."<<std::endl;
  }
  if (!partons->empty()) {
    partons->erase(partons->begin(),partons->end());
    msg.Debugging()<<"   Erased and deleted partons."<<std::endl;
  }
  me_ps_interface->PrepareTrees();
  msg.Debugging()<<"   Prepared Trees."<<std::endl;
  Flow::ResetCounter();
  msg.Debugging()<<"Out Apacic::CleanUpEvent() :"<<std::endl;
}









void Apacic::FillBlob(Blob * blob) {
  msg.Debugging()<<"##################################################"<<std::endl
		 <<"In Apacic::FillBlob() :"<<std::endl
		 <<"   "<<AP_proc->Nin()<<" -> "<<AP_proc->Nout()<<" process."<<std::endl;
    
  msg.Debugging()<<"Partons : "<<partons<<" "<<partons->size()<<std::endl;

  AMATOOLS::Vec4D cms = AMATOOLS::Vec4D(0.,0.,0.,0.);
  AMATOOLS::Vec4D pos = AMATOOLS::Vec4D(0.,0.,0.,0.);
    
  Parton * newp;
  for (int i=0;i<AP_proc->Nin();++i) {
    newp = new Parton(partons->size(),AP_proc->Flavs()[i],AP_proc->Momenta()[i]);
    newp->SetInfo('I');
    newp->SetNumber( partons->size() );
    newp->SetDec(blob);
    blob->AddToInPartons(newp);
    partons->push_back(newp);
    cms             = cms + AP_proc->Momenta()[i];
  }
  for (int i=AP_proc->Nin();i<AP_proc->Nin()+AP_proc->Nout();++i) {
    newp   = new Parton(partons->size(),AP_proc->Flavs()[i],AP_proc->Momenta()[i]);
    newp->SetInfo('H');
    newp->SetNumber( partons->size() );
    newp->SetProd(blob);
    partons->push_back(newp);
    blob->AddToOutPartons(newp);
  }
  
  blob->SetCMS(cms);
  blob->SetPosition(pos);
  blob->SetType(std::string("Hard Internal Partons (APACIC++ 2.0)"));
  
  msg.Debugging()<<"Out Apacic::FillBlob() :"<<std::endl;
  for (int i=0;i<AP_proc->Nin();++i) 
    msg.Debugging()<<blob->InParton(i)->Flav()<<" ";
  for (int i=AP_proc->Nin();i<AP_proc->Nin()+AP_proc->Nout();++i) 
    msg.Debugging()<<blob->OutParton(i-blob->NInP())->Flav()<<" ";
  msg.Debugging()<<std::endl;
  msg.Debugging()<<"########################################################"<<std::endl;
}

