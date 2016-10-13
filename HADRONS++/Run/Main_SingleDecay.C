#include "Main.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Decays/Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"

#ifdef USING__ROOT
static map<Hadron_Decay_Channel*, TFile*> rootfiles;
static map<Hadron_Decay_Channel*, vector<TH1D*> > TwoInvMassHists;
static map<Hadron_Decay_Channel*, vector<TH1D*> > ThreeInvMassHists;
static map<Hadron_Decay_Channel*, vector<TH1D*> > ThetaHists;
static map<Hadron_Decay_Channel*, vector<TH1D*> > EnergyHists;
#endif
static Flavour mother_flav;
static SHERPA::Hadron_Decay_Handler* hadrons;


void InitialiseGenerator(int argc, char *argv[])
{
  if(argc<2) {
    cout<<"Usage: ./SingleDecay <PDG_CODE>"<<endl;
    THROW(normal_exit,"you didn't specify the decaying particle by PDG code.");
  }

  small_sherpa_init(argc, argv);

  hadrons = new SHERPA::Hadron_Decay_Handler(".", "Fragmentation.dat");

  mother_flav = Flavour( (kf_code) abs(ToType<int>(argv[1])) );
  mother_flav.SetStable(false);
  if(ToType<int>(argv[1])<0) mother_flav=mother_flav.Bar();
  if(hadrons->DecayMap()->FindDecay(mother_flav)==NULL)
    THROW(fatal_error, "Didn't find "+ToString<Flavour>(mother_flav)+
	  " in HadronDecays.dat.");
  
  // set all decay channel BR's equal, such that we generate the same amount of
  // each decay channel to be tested
  PHASIC::Decay_Table* table=hadrons->DecayMap()->FindDecay(mother_flav);
  for(size_t i(0);i<table->size();++i)
    table->UpdateWidth(table->at(i), 1.0);
  //table->UpdateWidth(table->at(0), 1.0);
  //PRINT_VAR(table->at(0)->Name());
  
  rpa->gen.SetEcms(mother_flav.HadMass());
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;
}


void InitialiseAnalysis()
{
#ifdef USING__ROOT
  map<string, string> tags=Read_Write_Base::GlobalTags();
  string filepiece;
  map<string, string>::const_iterator it=tags.find("TAG_FILE_PIECE");
  if(it!=tags.end())
    filepiece="."+it->second;
  PHASIC::Decay_Table* table=hadrons->DecayMap()->FindDecay(mother_flav);
  for(int i(0);i<table->size();++i) {
    PHASIC::Decay_Channel* dc=table->at(i);
    Hadron_Decay_Channel* hdc = dynamic_cast<Hadron_Decay_Channel*>(dc);
    string fname = hdc->FileName();
    rootfiles[hdc]=
      new TFile(("Analysis/"+fname.erase(fname.length()-4)+filepiece+".root").c_str(),"RECREATE");

    TH1D* hist;
    int currenthist=0;
    int noutp(hdc->NOut());
    for(int i=0; i<noutp; i++) {
      Flavour flav = hdc->GetDecayProduct(i);
      string name = "costheta_"+flav.IDName()+"_"+ToString(i);
      string xtitle = "cos(#Theta) of "+flav.RootName();
      string ytitle = "#frac{1}{#Gamma} #frac{d#Gamma}{dcos(#Theta)}";
      hist = makeTH1D(name, "", 50, -1.0, 1.0, xtitle, ytitle);
      ThetaHists[hdc].push_back(hist);
      currenthist++;
    }
  
    currenthist=0;
    for(int i=0; i<noutp; i++) {
      Flavour flav = hdc->GetDecayProduct(i);
      string name = "E_"+flav.IDName()+"_"+ToString(i);
      string xtitle = "E of "+flav.RootName()+" [GeV]";
      string ytitle = "#frac{1}{#Gamma} #frac{d#Gamma}{dE} [GeV^{-1}]";
      double M = mother_flav.HadMass();
      double othermass = 0.0;
      for(int k=0; k<noutp; k++) {
	if(k!=i) othermass+=hdc->GetDecayProduct(k).HadMass();
      }
      double Emin = 0.9*flav.HadMass();
      double Emax = 1.1/2.0/M*(sqr(M)+sqr(flav.HadMass())-othermass);
      hist = makeTH1D(name, "", 50, Emin, Emax, xtitle, ytitle);
      EnergyHists[hdc].push_back(hist);
      currenthist++;
    }

    currenthist=0;
    for(int i=0; i<noutp; i++) {
      if(noutp<3) break;
      for(int j=i+1; j<noutp; j++) {
	Flavour flav1 = hdc->GetDecayProduct(i);
	Flavour flav2 = hdc->GetDecayProduct(j);
	string name = "q2_"+flav1.IDName()+"_"+flav2.IDName()+"_"+ToString(i)+ToString(j);
	string xtitle = "Inv. Mass q^{2}=(p_{"+flav1.RootName()+"}+p_{"
	  +flav2.RootName()+"})^{2} [GeV^{2}]";
	string ytitle = "#frac{1}{#Gamma} #frac{d#Gamma}{dq^{2}} [GeV^{-2}]";
	double othermass = 0.0;
	for(int k=0; k<noutp; k++) {
	  if(k==i || k==j) continue;
	  othermass+=hdc->GetDecayProduct(k).HadMass();
	}
	double q2min = 0.0;
	double q2max = 1.1*sqr(hdc->GetDecaying().HadMass()-othermass);
	hist = makeTH1D(name.c_str(), "", 50, q2min, q2max, xtitle, ytitle);
	TwoInvMassHists[hdc].push_back(hist);
	currenthist++;
      }
    }

    currenthist=0;
    for(int i=0; i<noutp; i++) {
      if(noutp<4) break;
      for(int j=i+1; j<noutp; j++) {
	for(int k=j+1; k<noutp; k++) {
	  Flavour flav1 = hdc->GetDecayProduct(i);
	  Flavour flav2 = hdc->GetDecayProduct(j);
	  Flavour flav3 = hdc->GetDecayProduct(k);
	  string name = "q2_"+flav1.IDName()+"_"+flav2.IDName()+"_"+flav3.IDName()
	    +"_"+ToString(i)+ToString(j)+ToString(k);
	  string xtitle = "Inv. Mass q^{2}=(p_{"+flav1.RootName()+
	    "}+p_{"+flav2.RootName()+"}+p_{"+flav3.RootName()+"})^{2} [GeV^{2}]";
	  string ytitle = "#frac{1}{#Gamma} #frac{d#Gamma}{dq^{2}} [GeV^{-2}]";
	  double othermass = 0.0;
	  for(int l=0; l<noutp; l++) {
	    if(l==i || l==j || l==k) continue;
	    othermass+=hdc->GetDecayProduct(k).HadMass();
	  }
	  double q2min = 0.0;
	  double q2max = 1.1*sqr(hdc->GetDecaying().HadMass());
	  hist = makeTH1D(name.c_str(), "", 50, q2min, q2max, xtitle, ytitle);
	  ThreeInvMassHists[hdc].push_back(hist);
	  currenthist++;
	}
      }
    }
  }
#endif
}


Blob_List* GenerateEvent()
{
  Blob_List* blobs = new Blob_List();
  Particle* mother_part = new Particle( 1,mother_flav,Vec4D(mother_flav.HadMass(),0.,0.,0.) );
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.HadMass());
  
  Blob* blob = blobs->AddBlob(btp::Hadron_Decay);
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->AddToInParticles(mother_part);

  try {
    hadrons->FillOnshellDecay(blob, NULL);
  } catch (Return_Value::code ret) {
    msg_Error()<<METHOD<<" Something went wrong for event: "<<*blobs
        <<endl;
    return blobs;
  }

  hadrons->CleanUp();

  msg_Events()<<*blobs<<std::endl;
  return blobs;
}


void AnalyseEvent(Blob_List* blobs)
{
  Blob* blob = blobs->FindFirst(btp::Hadron_Decay);
  if(blob==NULL) return;
#ifdef USING__ROOT
  Hadron_Decay_Channel* hdc=(*blob)["dc"]->Get<Hadron_Decay_Channel*>();

  int currenthist=0;
  for(int i=0; i<blob->NOutP(); i++) {
    double costheta = blob->OutParticle(i)->Momentum().CosTheta();
    ThetaHists[hdc][currenthist]->Fill(costheta);
    currenthist++;
  }
  
  currenthist=0;
  for(int i=0; i<blob->NOutP(); i++) {
    double energy = blob->OutParticle(i)->Momentum()[0];
    EnergyHists[hdc][currenthist]->Fill(energy);
    currenthist++;
  }

  currenthist=0;
  for(int i=0; i<blob->NOutP(); i++) {
    if(blob->NOutP()<3) break;
    for(int j=i+1; j<blob->NOutP(); j++) {
      Vec4D mom1 = blob->OutParticle(i)->Momentum();
      Vec4D mom2 = blob->OutParticle(j)->Momentum();
      double q2 = (mom1+mom2).Abs2();
      TwoInvMassHists[hdc][currenthist]->Fill(q2);
      currenthist++;
    }
  }

  currenthist=0;
  for(int i=0; i<blob->NOutP(); i++) {
    if(blob->NOutP()<4) break;
    for(int j=i+1; j<blob->NOutP(); j++) {
      for(int k=j+1; k<blob->NOutP(); k++) {
        Vec4D mom1 = blob->OutParticle(i)->Momentum();
        Vec4D mom2 = blob->OutParticle(j)->Momentum();
        Vec4D mom3 = blob->OutParticle(k)->Momentum();
        double q2 = (mom1+mom2+mom3).Abs2();
        ThreeInvMassHists[hdc][currenthist]->Fill(q2);
        currenthist++;
      }
    }
  }
#endif
}


void CleanUpEvent(Blob_List* blobs)
{
  blobs->Clear();
  Blob::Reset();
  Particle::Reset();
  Flow::ResetCounter();
  delete blobs;
}


void FinishGenerator()
{
  delete hadrons;
}


void FinishAnalysis()
{
#ifdef USING__ROOT
  map<Hadron_Decay_Channel*, TFile*>::iterator mit;
  for(mit=rootfiles.begin();mit!=rootfiles.end();++mit) mit->second->Write();
#endif
}
