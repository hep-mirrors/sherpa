#include "Sherpa.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"

#include "Process_Base.H"

#include "Flavour.H"
#include "Random.H"


#include "Amegic.H"

#include "Sample_Analysis.H"
#include "Shower_Observables.H"

#include <algorithm>

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#include "Random.H"
#endif

namespace SHERPA {
  Sherpa mm;
}

// interface for calling from fortran or c:
extern "C" {
  void apainit_() {
    SHERPA::mm.Init();
    SHERPA::mm.CrossSections();
    //    mm.GenerateEvents();
  }
  void aparun_() {
    SHERPA::mm.OneEvent();
  } 
  void apaend_() {
    SHERPA::mm.Finalize();
  } 
}


extern "C" {
  void fpyshower_(int &);
  void fhawface_(int&, int*, double*);
  void altest_(double&,double&,double&);
}

using namespace SHERPA;
using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

//const double facycut=.125;
//const double facycut=.125;

bool Sherpa::Init(int argc,char* argv[]) {
  ParticleInit("./");
  rpa.Init("./",argc, argv);   // use comand line parameter

  if (!as) as       = new Running_AlphaS();
  if (!aqed) aqed   = new Running_AlphaQED();




  double sy = rpa.test.FactorYcut()*rpa.integ.Ycut()*sqr(rpa.gen.Ecms());
  double sa = sqr(rpa.gen.Ecms());

  cout<<" a = "<<(*as)(sy)<<endl;
  cout<<" b = "<<(*as)(sa)<<endl;
  //  exit(0);

  tune    = 0;  // *AS* default is 0


  partons           = new Parton_List;
  blobs             = new Blob_List;
  AM                = new Amegic("./"); // ,isr;

  int jetnumber     = 5;   // *AS* WARNING this influences the jetveto! Have to apply
  // modified sudakovs as well!
  
  isr=0; // NO Initial-Shower (isr has to be initialised by Amegic)
  hard_interface    = new ME_PS_Interface(isr,jetnumber); 
  cout<<" init Soft Interface "<<endl;
  soft_interface    = new Soft_Interface();
  cout<<" init Soft Interface done"<<endl;

  return 1;
}

Sherpa::~Sherpa() {
  if (blobs) CleanUpEvent();

  if (hard_interface) delete hard_interface;
  if (soft_interface) delete soft_interface;
  if (AM)             delete AM;
  if (as)             delete as;     as=0;    // global vars
  if (aqed)           delete aqed;   aqed=0;  // global vars
}


struct Rescale_Data {
  int    njet;
  double xstot;
  double r;
  double r_dir;
  double r_res1;
  double r_res2;
  double r_nat;
  Process_Base * proc;
};

class Compare_Rescale_Data {
public:
  int operator()(const Rescale_Data & a, const Rescale_Data & b) {
    if (a.njet<b.njet) return 1;
    return 0;
  }
};


bool Sherpa::RescaleJetrates() {
  Process_Base * procs = AM->Processes();
  // we have somehow to determine wich process is included in which other process!
  // ???!!!
  // assuming just 2,3,4 jet events

  // fill list of totalxs and number of outgoing particles
  vector<Rescale_Data> pd(procs->Size());
  for (int i=0; i<procs->Size();++i) {
    pd[i].xstot = (*procs)[i]->Total();
    pd[i].njet  = (*procs)[i]->Nout();
    pd[i].proc  = (*procs)[i];
  }
  int in2   = -1;
  for (int i=0; i<pd.size();++i) {
    if (pd[i].njet==2) {
      if (in2!=-1) {
	msg.Out()<<" WARNING: someting wrong in RescaleJetrates: two 2->2 process groups!"<<endl;
      }
      in2=i;
    }
  }  

  if (in2!=-1) {
    if (pd.size()>=2) {
      // sort
      std::sort(pd.begin(),pd.end(),Compare_Rescale_Data());

      // check that njet==n is followed by njet==n+1
      int ok=1;
      for (int i=1; i<pd.size();++i) {
	if (pd[i-1].njet!=pd[i].njet-1) {
	  msg.Out()<<" someting wrong in RescaleJetrates: not a correct list of processes"<<endl;
	}
      }

      // determine r=sigma_n/sigma_2
      double sum_res2=0.;
      double sum_res1=0.;
      double sum_nat =0.;
      double sum_dir =0.;
      double sum_r   =0.;
      for (int i=pd.size()-1; i>=0;--i) {
	pd[i].r_nat=pd[i].r_res1=pd[i].r_dir=pd[i].r=pd[i].xstot/pd[0].xstot;
	if (i!=pd.size()-1) {
	  pd[i].r_nat -=pd[i+1].r;
	  pd[i].r_res1-=sum_dir;
	  if (i==0) {
	    pd[0].r_dir  = 1. - sum_dir;
	    pd[0].r_res1 = 1. - sum_res1;
	  }
	}
	sum_r       +=pd[i].r;                    // sig(i)/(sum sig(j)) cf. (**)
	sum_dir     +=pd[i].r_dir;                // sig(i)/sig(1) i!=1  
	sum_nat     +=pd[i].r_nat;                // [sig(i)-sig(i+1)]/sig(1)
	sum_res1    +=pd[i].r_res1;               // res1
	sum_res2    +=pd[i].r_res2;               // not implemented
      }

      // rescale and output
      for (int i=0; i<pd.size(); ++i) {
	//	pd[i].proc->RescaleXSec(pd[i].r_nat/pd[i].r);
      }

      for (int i=0;i<pd.size();++i) { 
	pd[i].r*=1./sum_r;                        //  (**)

	msg.Out()<<" "<<pd[i].njet<<"-Jet  sigma="<<pd[i].xstot<<endl;
	msg.Out()<<"\t\t"<<pd[i].r<<"\t"<<pd[i].r_nat<<"\t"<<pd[i].r_dir<<"\t"
		 <<pd[i].r_res1<<"\t"<<pd[i].r_res2<<endl;
	
      }
      // actualize "All_Processes"
      double dummy=1.;
      procs->RescaleXSec(dummy);
      
      msg.Out()<<" Sum: "<<endl;
      msg.Out()<<"\t\t"<<sum_r<<"\t"<<sum_nat<<"\t"<<sum_dir<<"\t"
	       <<sum_res1<<"\t"<<sum_res2<<endl;
    }
  }
  else {
    msg.Out()<<" someting wrong in RescaleJetrates: no 2->2 process group! "<<endl;
    msg.Out()<<"  no rescaling !"<<endl;
  }

}

bool Sherpa::CrossSections() {
  msg.SetLevel(0);

  if (AM->InitializeProcesses()) { 
    msg.Out()<<" Matching: yf   = "<<rpa.test.FactorYcut()<<endl;
    msg.Out()<<" Matching: nllf = "<<rpa.test.FactorNLLQ()<<endl;


    if (tune) { 
      int  mode_dir = 448;
      mkdir("Tuning",mode_dir); 
      AM->SetResDir(std::string("Tuning"));
      // set scale for alphaS factors
      if (rpa.me.KFactorScheme()==1 && rpa.me.ScaleScheme()==2)
	AM->Processes()->SetScale(rpa.test.FactorYcut()*rpa.integ.Ycut()*sqr(rpa.gen.Ecms()));
      if (!AM->LookUpXSec(rpa.integ.Ycut(),1,std::string("dY_cut"))) return 0;
      return 1;
    }
    else {
      // set scale for alphaS factors
      if (rpa.me.KFactorScheme()==1 && rpa.me.ScaleScheme()==2)
	AM->Processes()->SetScale(rpa.test.FactorYcut()*rpa.integ.Ycut()*sqr(rpa.gen.Ecms()));
      if (AM->CalculateTotalXSec()) {
	// determine (rescaled) jetrates
	RescaleJetrates();
	// successful
	return 1;
      }
    }
  }
  return 0;
}

bool Sherpa::GenerateEvents() {
  if (rpa.gen.NumberOfEvents()==0) return 1;
  if (MPIGenerateEvents()) return 1;
  msg.SetLevel(rpa.gen.Output());
  //  RunPythiaTest(); 

  msg.Tracking()<<"Start generating events : "<<std::endl
		<<"Use Amegic for the hard processes in event generation."<<std::endl;

  Process_Base * proc;
  Blob         * blob;

  //---- init a sample analyse 
  Sample_Analysis analysis;
  analysis.Init();

  //  Ran.ReadInStatus("RandomA.dat",169);

  msg.Out()<<" Starting event generation now. "<<std::endl;
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) {
    if (n%2500==0) {
      //      int no=Ran.WriteOutStatus("RandomC.dat");
      int no=0;
      msg.Out()<<" "<<n<<" Events ( ran="<<no<<")"<<std::endl;
//       cout<<" Press Enter to continue ! "<<endl;
//       string enter;
//       cin>>enter;
    }
    msg.Events()<<"##########################################################"<<std::endl
		<<"#############################"<<n<<"th event #############"<<std::endl
		<<"##########################################################"<<std::endl;

    // determine one hard event, with kinematics (and an factor [as(s*y)/as(s)]^(nout-2) )
    if (AM->Processes()->OneEvent()) {
      proc = AM->Processes()->Selected();
      int stat=1;
      do { 
	if (stat==3) {
	  // need new kinematics
	  AM->Processes()->SameEvent();
	}

	blob = new APHYTOOLS::Blob();
	FillBlob(blob,proc);
	blob->SetId(blobs->size());
	blobs->push_back(blob);

	blob->BoostInCMS();

	// determine sudakov factors	
	if (!(hard_interface->Treat(proc,blob,1))) {
	  // event rejected due to Sudakov-n-AlphaS-weight
	  stat=4;
	}
	else {
	  // perform parton shower
	  stat=hard_interface->PerformShower(proc,1);
	  if ((stat==3)  ) {
   	    // rejection of kinematics due to minjet-veto
	    //	    cout<<" stat==3"<<endl;
	    CleanUpEvent();
	  }
	}
      } while (stat==3);

     analysis.AfterME(blobs);

     if ((stat==4)  ) {
       // event rejected due to Sudakov-n-AlphaS-weight	
       --n;
       CleanUpEvent();
       continue;
     }

      if ((stat==0)  ) {
	msg.Error()<<"ERROR in Sherpa::GenerateEvents"<<std::endl
		   <<"   The hard interface did not know how to perform the shower."
		   <<std::endl;
	--n;
	CleanUpEvent();
	continue;
      }


      hard_interface->ExtractPartons(blobs,partons);
      msg.SetPrecision(4);
      msg.Events()<<(*blobs)<<std::endl;
      msg.SetPrecision(12);

      analysis.AfterPartonShower(blobs);

      soft_interface->PerformFragmentation(blobs,partons);
//          analysis.AfterHadronization(blobs);

      msg.SetPrecision(4);
      msg.Events()<<(*blobs)<<std::endl;
      msg.SetPrecision(12);


      CleanUpEvent();
    }
    else {
      msg.Error()<<" Warning: AM->Processes()->OneEvent() in Sherpa.C failed "<<std::endl;
    }
  }

  // Analysis: write out histos
  analysis.Finish();

  return 1;
}


void Sherpa::OneEvent() {
  int stat=1;
  do {
    // remove last event
    CleanUpEvent();

    // determine one hard event, with kinematics (and an factor [as(s*y)/as(s)]^(nout-2) )
    if (AM->Processes()->OneEvent()) {
      Process_Base * proc = AM->Processes()->Selected();

      stat=1;
      do { 
	if (stat==3) {
	  // need new kinematics
	  AM->Processes()->SameEvent();
	}

	Blob * blob= new Blob();
	FillBlob(blob,proc);

	blob->SetId(blobs->size());
	blobs->push_back(blob);
	blob->BoostInCMS();

	// determine sudakov factors	
	if (!(hard_interface->Treat(proc,blob,1))) {
	  // event rejected due to Sudakov-n-AlphaS-weight
	  stat=4;
	  //	  continue; 
	}
	else {
	  // perform parton shower
	  stat=hard_interface->PerformShower(proc,1);
	  if ((stat==3)  ) {
   	    // rejection of kinematics due to minjet-veto
	    //	    cout<<" stat==3"<<endl;
	    CleanUpEvent();
	  }
	}
      } while (stat==3);

      if ((stat==4)  ) {
	// event rejected due to Sudakov-n-AlphaS-weight	
	continue;
      }

      if ((stat==0)  ) {
	msg.Error()<<"ERROR in Sherpa::GenerateEvents"<<std::endl
		   <<"   The hard interface did not know how to perform the shower."
		   <<std::endl;
	continue;
      }

      hard_interface->ExtractPartons(blobs,partons);

      soft_interface->PerformFragmentation(blobs,partons);

      msg.SetPrecision(4);
      msg.Events()<<(*blobs)<<std::endl;
      msg.SetPrecision(12);

    }
    else {
      msg.Error()<<" Warning: AM->Processes()->OneEvent() in Sherpa.C failed "<<std::endl;
    }
  } while (stat!=1);
}

void Sherpa::Finalize() {
  // nothing to be done
}


void Sherpa::RunPythiaTest() {
  bool calc_alphas=0;
  int process=1;       // 0 - ee-> jets;  1 - drell-yan

  if (calc_alphas) {
    msg.Out()<<"calculating alphas with Pythia: "<<std::endl;
    double q2_max=sqr(91.2);
    double q2_min=sqr(.912);
    int    n =10;
    for (int i=0;i<=n;++i) {
      double q2=q2_min*pow((q2_max/q2_min),double(i)/double(n));
      double alphas=0, lambda=0;
      altest_(q2,lambda,alphas);
      cout<<" "<<q2<<" \t"<<alphas<<" \t"<<sqr(lambda)<<endl;
    }
    //  exit(0);
  }
  

  msg.Out()<<"Start generating events with Pythia: "<<std::endl;


  const int maxentries = 2000;

  double  * pjet = new double[4*maxentries];
  int       nk;
  int * kfjet  = new int[maxentries];

  Parton_List pl;


  Primitive_Analysis ana;
  ana.AddObservable(new Jetrates(11,1.e-6,1.,180,0));
  ana.AddObservable(new Multiplicity(00,-0.5,50.5,51,0));




  double E  =rpa.gen.Ecms()*0.5;


  msg.Out()<<" Starting event generation now. "<<std::endl;
  for (int n=1;n<=rpa.gen.NumberOfEvents();n++) {
    if (n%2500==0) {
      msg.Out()<<"event ="<<n<<endl;
    }
    //    msg.Out()<<"event ="<<n<<endl;

    // call pythia
    fpyshower_(n);

    //    msg.Out()<<"event done ="<<n<<endl;

    // translate to pl
    
    fhawface_(nk,kfjet,pjet);
    
    //    msg.Out()<<"nk ="<<nk<<endl;

    if (process==0) {
      pl.push_back(new Parton(0,Flavour(kf::e),Vec4D(E,0,0,E)));
      pl.push_back(new Parton(1,Flavour(kf::e).Bar(),Vec4D(E,0,0,-E)));
    }
    else if (process==1) {
      // pseudo initiator!
      pl.push_back(new Parton(0,Flavour(kf::gluon),Vec4D(E,0,0,E)));
      pl.push_back(new Parton(1,Flavour(kf::gluon),Vec4D(E,0,0,-E)));
    }
    for (int i=0; i<nk; ++i) {
      Flavour flav = Flavour(kf::code(abs(*(kfjet+i))));
      if ((*(kfjet+i))<0) flav = flav.Bar();
      Vec4D mom;
      for(int j=0; j<4; ++j) mom[j] = *(pjet+i+j*2000);

      if (flav.IsLepton()) {
	msg.Events()<<" ignoring "<<flav<<" "<<mom<<endl;
      }
      else if (flav.IsDiQuark() || flav.IsHadron()) {
	msg.Events()<<" assuming Beam remnant "<<flav<<" "<<mom<<endl;
	if (mom[3]>0) {
	  mom=pl[0]->Momentum()-mom;
	  pl[0]->SetMomentum(mom);
	}
	else {
	  mom=pl[1]->Momentum()-mom;
	  pl[1]->SetMomentum(mom);
	}
      } 
      else {
	Parton * parton = new Parton(pl.size(),flav,mom);
	parton->SetStatus(1);
	
	pl.push_back(parton);
      }
    }

    msg.Events()<<" pln="<<pl.size()<<endl;
    for (int i=0; i<pl.size();++i) { 
      msg.Events()<<i<<" :"<<pl[i]<<endl;
    }

    
    // do analysis
    ana.DoAnalysis(pl,1.);

    // delete partons
    for (int i=0; i<pl.size();++i) {
      delete pl[i];
    }
    pl.clear();
  }

  // Analysis: write out histos
  ana.FinishAnalysis("pythia_HE",0);

    
  delete kfjet;
  delete pjet;

  exit(0);

}


void Sherpa::CleanUpEvent() {
  if (!blobs->empty()) {
    // delete Blobs
    for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) delete (*blit);
    //blobs->erase(blobs->begin(),blobs->end());
    blobs->clear();
  }

  if (!partons->empty()) partons->erase(partons->begin(),partons->end());
  //count->ResetCounter();
  Flow::ResetCounter();
}

void Sherpa::FillBlob(Blob * blob,Process_Base * proc) {
  msg.Debugging()<<"##################################################"<<std::endl
		 <<"In Sherpa::FillBlob() :"<<std::endl
		 <<"   "<<proc->Nin()<<" -> "<<proc->Nout()<<" process."<<std::endl;
    
  msg.Debugging()<<"Partons : "<<partons<<" "<<partons->size()<<std::endl;

  AMATOOLS::Vec4D cms = AMATOOLS::Vec4D(0.,0.,0.,0.);
  AMATOOLS::Vec4D pos = AMATOOLS::Vec4D(0.,0.,0.,0.);
    
  Parton * newp;
  for (int i=0;i<proc->Nin();++i) {
    newp = new Parton(partons->size(),proc->Flavs()[i],proc->Momenta()[i]);
    newp->SetInfo('I');
    newp->SetNumber( partons->size() );
    newp->SetDec(blob);
    blob->AddToInPartons(newp);
    partons->push_back(newp);
    cms             = cms + proc->Momenta()[i];
    //    cout<<" AddI :"<<newp<<endl;
  }
  for (int i=proc->Nin();i<proc->Nin()+proc->Nout();++i) {
    newp   = new Parton(partons->size(),proc->Flavs()[i],proc->Momenta()[i]);
    newp->SetInfo('H');
    newp->SetNumber( partons->size() );
    newp->SetProd(blob);
    partons->push_back(newp);
    blob->AddToOutPartons(newp);
    //    cout<<" AddO :"<<newp<<endl;
  }
  blob->SetCMS(cms);
  blob->SetPosition(pos);
  blob->SetType(std::string("Hard ME (AMEGIC++2.0)"));
  
  msg.Debugging()<<"Out Sherpa::FillBlob() :"<<std::endl;
  for (int i=0;i<proc->Nin();++i) 
    msg.Debugging()<<blob->InParton(i)->Flav()<<" ";
  for (int i=proc->Nin();i<proc->Nin()+proc->Nout();++i) 
    msg.Debugging()<<blob->OutParton(i-blob->NInP())->Flav()<<" ";
  msg.Debugging()<<std::endl;
  msg.Debugging()<<"########################################################"<<std::endl;
}

bool Sherpa::MPIInit() {
#ifdef _USE_MPI_

  return 1;
#else
  return 0;
#endif

}




bool Sherpa::MPIGenerateEvents() {
#ifdef _USE_MPI_
  int mpi_rank = MPI::COMM_WORLD.Get_Rank();
  int mpi_size = MPI::COMM_WORLD.Get_size();

  if (mpi_size<2) return 0;

  Ran.Set_Seed(-3784-mpi_rank*7);

  // --- define MPI Messages --- 
  MPI_Parton    sample_parton;
  MPI::Datatype mpi_parton_type;
  MPI::Datatype type[5]     = {MPI::INT,MPI::CHAR,MPI::INT,MPI::DOUBLE,MPI::INT};
  int           blocklen[5] = {2       ,1        ,3       ,4          ,2};
  MPI::Aint     disp[5];

  disp[0] = MPI::Get_address(&sample_parton.id);
  disp[1] = MPI::Get_address(&sample_parton.m_info);
  disp[2] = MPI::Get_address(&sample_parton.m_fl);
  disp[3] = MPI::Get_address(&sample_parton.m_mom);
  disp[4] = MPI::Get_address(&sample_parton.blob_owned_by);
  for (short int i=4;i>=0;i--) disp[i] -= disp[0];

  mpi_parton_type = MPI::Datatype::Create_struct(5,blocklen,disp,type);
  mpi_parton_type.Commit();
  // --- end define Message ---

  msg.Out()<<" Starting node "<<(mpi_rank+1)<<"/"<<mpi_size<<endl;
  //  if (mpi_rank==0) MPIInit();


  msg.Tracking()<<"Start generating events : "<<std::endl
		<<"Use Amegic for the hard processes in event generation."<<std::endl;

  Process_Base * proc;
  Blob         * blob;


  MPI::Status    status;
  const int NMAX=500;
  MPI_Parton buffer[NMAX];
  
  int tag_pl = 1;
  int tag_fin = 9;
  int dest_master = 0;


  // MASTER
  if (mpi_rank==0) {

    //---- init a test  analyse 
    Primitive_Analysis ana;
    ana.AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));
    msg.Out()<<" Starting event generation now. "<<std::endl;


    for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) {
      if (n%2500==0) 
	msg.Out()<<" "<<n<<" Events "<<std::endl;


      // collect information of a SLAVE !!!
      while (!MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG)) {
	// wait
	double a=1.3;
	double b=1.4;
	a*b;
      }

      MPI::COMM_WORLD.Recv(buffer,NMAX,mpi_parton_type, MPI::ANY_SOURCE,MPI::ANY_TAG , status);
      int np = status.Get_count(mpi_parton_type);
      if (np!=MPI::UNDEFINED) {
	//	cout<<"master:"<<endl<<" np="<<np<<endl;
	for (int j=0;j<np;++j) {
	  if (buffer[j].blob_owned_by>=blobs->size()) {
	    // generate new blob
	    blob = new APHYTOOLS::Blob();
	    blob->SetId(blobs->size());
	    if (buffer[j].blob_owned_by==0)  blob->SetType(std::string("Hard ME (AMEGIC++2.0)")); 
	    else                             blob->SetType(std::string("FS Parton Shower (APACIC++2.0)"));
	    blobs->push_back(blob);
	  }
	  Parton * p = MPI2Parton(buffer[j]);
	  if (buffer[j].is_initial)    blob->AddToInPartons(p);
	  else                         blob->AddToOutPartons(p);
	}
      }
      else {
	cout<<" error in transmission "<<endl;
      }

      //      cout<<(*blobs)<<endl;


      ana.DoAnalysis(*blobs,1.);
      CleanUpEvent();

    }

    // send break signals !!!!
    for (int i=1;i<mpi_size;++i) {
      int fin=1;
      MPI::COMM_WORLD.Send(&fin,1,MPI::INT,i,tag_fin); 
    }


    // empty mpi buffers
    while (MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG)) {
      MPI::COMM_WORLD.Recv(buffer,NMAX,mpi_parton_type, MPI::ANY_SOURCE,MPI::ANY_TAG , status);
      cout<<" destroy unwanted data from SLAVE "<<status.Get_source()<<endl;
    }

    // Analysis: write out histos
    ana.FinishAnalysis("testout4jet_c",0);


  }
  else {
    // SLAVE (mpi_rank!=0)
    int break_signal=0;

    long int nevt=rpa.gen.NumberOfEvents()/(mpi_size-1)+1;
    long int no  =0;

    for (;;) {

      // determine one hard event, with kinematics (and an factor [as(s*y)/as(s)]^(nout-2) )
      if (AM->Processes()->OneEvent()) {
	proc = AM->Processes()->Selected();
	int stat=1;
	//	cout<<" #"<<no<<endl;
	do { 
	  if (stat==3) {
	    // need new kinematics
	    //	    cout<<" new kin"<<endl;
	    AM->Processes()->SameEvent();
	    //	    cout<<" done "<<endl;
	  }

	  blob = new APHYTOOLS::Blob();
	  FillBlob(blob,proc);
	  blob->SetId(blobs->size());
	  blobs->push_back(blob);
	  blob->BoostInCMS();
      
	  if (!(hard_interface->Treat(proc,blob,1))) {
	    // event rejected due to Sudakov-n-AlphaS-weight
	    stat=4;
	    //	    continue; 
	  }
	  else {
// 	  if (stat==3) {
// 	    cout<<" new blob "<<(*blob);
// 	  }

	  stat=hard_interface->PerformShower(proc,1);
	  if ((stat==3)  ) {
//  	    cout<<"NEW Feature in Sherpa::GenerateEvents"<<std::endl
//  		<<" got return value 3! dice again  "<<std::endl;
	 
// 	    cout<<" old blob"<<(*blob);
	    // *AS*  rejection test
	    // --n;
	    CleanUpEvent();
	  }
	  }
	} while (stat==3);

	if ((stat==4)  ) {
	  //	  cout<<" event rejected due to Sudakov-n-AlphaS-weight "<<endl;
	  //	  --n;
	  CleanUpEvent();
	  continue;
	}


	if ((stat==0)  ) {
	  msg.Error()<<"ERROR in Sherpa::GenerateEvents"<<std::endl
		     <<"   The hard interface did not know how to perform the shower."<<std::endl;
	  
	  // *AS*  rejection test
	  // --n;
	  CleanUpEvent();
	  continue;
	}

	hard_interface->ExtractPartons(blobs,partons);
	// *AS*      soft_interface->PerformFragmentation(blobs,partons);



	// send information to MASTER
	int j=0;
	for (Blob_Const_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
	  int nout=(*blit)->NOutP();
	  int nin=(*blit)->NInP();
	  for (int i=0;i<nin;++i) {
	    if (j<NMAX)
	      Parton2MPI((*blit)->InParton(i),buffer[j]);
	    buffer[j].is_initial=1;
	    buffer[j].blob_owned_by=(*blit)->Id();
	    ++j;
	  }
	  for (int i=0;i<nout;++i) {
	    if (j<NMAX)
	      Parton2MPI((*blit)->OutParton(i),buffer[j]);
	    buffer[j].is_initial=0;
	    buffer[j].blob_owned_by=(*blit)->Id();
	    ++j;
	  }
	}

	// check for BREAK signal
	if (MPI::COMM_WORLD.Iprobe(0, tag_fin)) {
	  int fin;
	  MPI::COMM_WORLD.Recv(&fin, 1, MPI::INT, 0, tag_fin);
	  cout<<"SLAVE "<<mpi_rank<<" received finish-tag after "<<no<<" events"<<endl;
	  break_signal=1;
	}
	else {
	  //	  cout<<" SLAVE "<<mpi_rank<<" sending "<<j<<" partons of bloblist "<<endl;
	  //	cout<<(*blobs)<<endl;
	  MPI::COMM_WORLD.Send(buffer,j,mpi_parton_type,dest_master,tag_pl);
	  ++no;
	}
	CleanUpEvent();

      }
      else {
	msg.Error()<<" Warning: AM->Processes()->OneEvent() in Sherpa.C failed "<<std::endl;
	
      }

      if (no>nevt) { 
	break_signal=1;
	cout<<" SLAVE "<<mpi_rank<<" has finished its share"<<endl;
      }

      if (break_signal) break;
    }
  }

  return 1;

#else
  return 0;
#endif
}

