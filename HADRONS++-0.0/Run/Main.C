#include "Exception.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Particle.H"
#include "Hadrons.H"
#include "Model_Base.H"
#include "Standard_Model.H"

#ifdef USING__ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h"
#endif

#include "Data_Reader.H"

using namespace std;
using namespace HADRONS;
using namespace ATOOLS;
using namespace MODEL;

int main(int argc, char *argv[])
{
#ifdef USING__ROOT
  if(argc!=2) {
    cout<<"Usage: ./Hadrons <PDG_CODE>"<<endl;
    abort();
  }
  set_terminate(Exception_Handler::Terminate);
  set_unexpected(Exception_Handler::Terminate);
  signal(SIGSEGV,Exception_Handler::SignalHandler);
  signal(SIGINT,Exception_Handler::SignalHandler);
  signal(SIGBUS,Exception_Handler::SignalHandler);
  signal(SIGFPE,Exception_Handler::SignalHandler);
  signal(SIGABRT,Exception_Handler::SignalHandler);
  signal(SIGTERM,Exception_Handler::SignalHandler);
  signal(SIGXCPU,Exception_Handler::SignalHandler);
  try {
    ATOOLS::ParticleInit("./");
    rpa.Init("./","Run.dat",argc,argv);
    Model_Base* model = new Standard_Model("./","Model.dat");
    rpa.gen.SetModel(model);
//   rpa.gen.SetEcms(91.2);
    msg.SetModifiable(true);
    Flavour mother_flav(kf::code(ToType<int>(argv[1])));
    msg.Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;
    Hadrons* hadrons = new Hadrons( "../../SHERPA-1.0/Run/Decaydata/",
                                    string("HadronDecays.dat"),
                                    string("HadronConstants.dat") ) ;

//     Hadrons* hadrons = new Hadrons( "./Decaydata/",
//                                     string("HadronDecays.dat"),
//                                     string("HadronConstants.dat") ) ;
    TApplication myapp(string("Hadrons").c_str(),&argc,&(*argv));

    if( !hadrons->FindDecay(mother_flav) || mother_flav.IsStable()) {
      msg.Error()<<"  The selected particle is either set to stable in Hadron.dat or "
        <<"it does not have a decay table in Decaydata."<<endl;
      abort();
    }
    int nevents;
    if(ATOOLS::rpa.gen.NumberOfEvents()==0) {
      cout<<"Enter the number of events: "<<endl;
      cin>>nevents;
    }
    else nevents=ATOOLS::rpa.gen.NumberOfEvents();

    TFile* rootfile; TH1D* hist;

    for(int i=0; i<nevents; i++) {

      Particle* mother_part = new Particle( 1,mother_flav,Vec4D(mother_flav.PSMass(),0.,0.,0.) );

      Blob* fragmentation_blob = new Blob();
      fragmentation_blob->AddToOutParticles(mother_part);
      fragmentation_blob->SetType(btp::Fragmentation);

      Blob_List* blob_list(new Blob_List());
      blob_list->push_back(fragmentation_blob);

      Blob* motherdecay_blob=hadrons->CreateDecayBlobSkeleton(mother_part,blob_list,NULL);
      hadrons->PerformDecay(motherdecay_blob,blob_list,NULL,NULL);
      msg.Events()<<"-------------------------------------------"<<endl;
      msg.Events()<<(*blob_list)<<endl;
      msg.Events()<<"-------------------------------------------"<<endl;

      if(i==0) {
        Hadron_Decay_Channel* hdc = (*motherdecay_blob)["hdc"]->Get<Hadron_Decay_Channel*>();
        string name     = hdc->ChannelName();
        const ATOOLS::Flavour* flavs = hdc->Flavours();
        string filename = flavs[0].ShellName()+"_";
        for(unsigned int j=1;j<hdc->NOut()+1;j++) {
          filename+=flavs[j].ShellName();
        }
        cout<<"Filling histogram for process "<<name<<endl;
        rootfile = new TFile(("Analysis/"+filename+".root").c_str(),"RECREATE");
        hist     = new TH1D("q2", ("q^{2} in "+name).c_str() ,100,0.0,4.0);
      }

      double q2 = (
                    motherdecay_blob->OutParticle(0)->Momentum()
                    + motherdecay_blob->OutParticle(1)->Momentum()
//                     + motherdecay_blob->OutParticle(3)->Momentum()
                  ).Abs2();
      hist->Fill(q2);

      fragmentation_blob->Reset();
      blob_list->Clear();
      delete blob_list;
      if(i%1000==0) msg.Info()<<i<<" events passed."<<std::endl;
    }
    msg.Info()<<"Generated "<<nevents<<" events"<<endl;
    delete hadrons;
    delete model;
    hist->Draw();
    myapp.Run(kTRUE);
    hist->Write();
    rootfile->Close();
  }
  catch (Exception Sherpaexception) {
    Sherpaexception.UpdateLogFile();
    msg.Error()<<Sherpaexception<<endl;
    terminate();
  }
  catch (std::exception stdexception) {
    cout<<"Sherpa: throws exception "
	     <<stdexception.what()<<" ..."<<endl;
    terminate();
  }
#endif
}
