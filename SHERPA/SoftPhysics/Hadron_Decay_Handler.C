#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include <algorithm>

using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;
using namespace HADRONS;
using namespace METOOLS;

Hadron_Decay_Handler::Hadron_Decay_Handler(string path, string fragfile) :
  Decay_Handler_Base()
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(path);
  dr.SetInputFile(fragfile);

  m_qedmode=dr.GetValue<size_t>("HADRON_DECAYS_QED_CORRECTIONS",1);
  double max_propertime = dr.GetValue<double>("MAX_PROPER_LIFETIME",-1.0);
  if( max_propertime > 0.0) {
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
	kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      if (flav.IsOn() && flav.IsHadron() && !flav.IsStable() &&
          0.197e-12>max_propertime*flav.Width() && flav.Kfcode()!=kf_K)
      {
        flav.SetStable(true);
      }
    }
  }

  string decaypath=dr.GetValue<string>("DECAYPATH",string("Decaydata/"));
  string decayfile=dr.GetValue<string>("DECAYFILE",string("HadronDecays.dat"));
  string decayconstfile=dr.GetValue<string>("DECAYCONSTFILE",
                                            string("HadronConstants.dat"));
  string bdecayfile=dr.GetValue<string>("B_DECAYFILE",
					string("Partonic_b/Decays.dat"));
  string cdecayfile=dr.GetValue<string>("C_DECAYFILE",
					string("Partonic_c/Decays.dat"));
  string aliasfile=dr.GetValue<string>("HADRONALIASESFILE",
                                       string("HadronAliases.dat"));
  string aliasdecayfile=dr.GetValue<string>("ALIASDECAYFILE",
                                            string("AliasDecays.dat"));
  m_mass_smearing=dr.GetValue<int>("SOFT_MASS_SMEARING",1);
  m_spincorr=rpa->gen.SoftSC();

  decaypath=rpa->gen.Variable("SHERPA_SHARE_PATH")+"/"+decaypath;
  My_In_File::OpenDB(decaypath);
  Hadron_Decay_Map * dmap = new Hadron_Decay_Map(this);
  dmap->ReadInConstants(decaypath, decayconstfile);
  dmap->ReadInPartonicDecays(Flavour(kf_b),decaypath,bdecayfile);
  dmap->ReadInPartonicDecays(Flavour(kf_c),decaypath,cdecayfile);
  dmap->ReadHadronAliases(decaypath, aliasfile);
  dmap->Read(decaypath, decayfile, true);
  dmap->Read(decaypath, aliasdecayfile);
  dmap->Initialise();
  dmap->ReadFixedTables(decaypath, "FixedDecays.dat");
  p_decaymap=dmap;
  
  p_mixinghandler = new Mixing_Handler();
  p_mixinghandler->SetModel(dmap->StartModel());
  dmap->SetMixingHandler(p_mixinghandler);
  My_In_File::CloseDB(decaypath);
}

Hadron_Decay_Handler::~Hadron_Decay_Handler()
{
  Hadron_Decay_Map* dmap=dynamic_cast<Hadron_Decay_Map*>(p_decaymap);
  delete dmap; p_decaymap=NULL;
  delete p_mixinghandler; p_mixinghandler=NULL;
}

void Hadron_Decay_Handler::TreatInitialBlob(ATOOLS::Blob* blob,
                                            METOOLS::Amplitude2_Tensor* amps,
                                            const Particle_Vector& origparts)
{
  if (blob->Has(blob_status::needs_showers)) return;
  if (RejectExclusiveChannelsFromFragmentation(blob)) return;
  Decay_Handler_Base::TreatInitialBlob(blob, amps, origparts);
}

Decay_Matrix* Hadron_Decay_Handler::FillDecayTree(Blob * blob, Spin_Density* s0)
{
  Blob* mixingblob=p_mixinghandler->PerformMixing(blob->InParticle(0));
  if (mixingblob) {
    Blob_List::iterator bit;
    for (bit=p_bloblist->begin();bit!=p_bloblist->end();++bit) 
      if (*bit==blob) p_bloblist->erase(bit);
    p_bloblist->push_back(mixingblob);
    CreateDecayBlob(mixingblob->OutParticle(0));
    return Decay_Handler_Base::FillDecayTree
      (mixingblob->OutParticle(0)->DecayBlob(), s0);
  }
  else return Decay_Handler_Base::FillDecayTree(blob, s0);
}

Amplitude2_Tensor*
Hadron_Decay_Handler::FillOnshellDecay(Blob *blob, Spin_Density* sigma)
{
  Amplitude2_Tensor* amps=Decay_Handler_Base::FillOnshellDecay(blob,sigma);
  Decay_Channel* dc=(*blob)["dc"]->Get<Decay_Channel*>();
  Hadron_Decay_Channel* hdc=dynamic_cast<Hadron_Decay_Channel*>(dc);
  //if (blob->InParticle(0)->Flav().IsB_Hadron())
  //  msg_Out()<<METHOD<<":\n"<<(*blob)<<"\n";
  if(hdc && !hdc->SetColorFlow(blob)) {
    msg_Error()<<METHOD<<" failed to set color flow, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  return amps;
}

void Hadron_Decay_Handler::CreateDecayBlob(Particle* inpart)
{
  DEBUG_FUNC(inpart->RefFlav());
  if(inpart->DecayBlob()) abort();
  if(!Decays(inpart->Flav())) return;
  if(inpart->Time()==0.0) inpart->SetTime();
  Blob* blob = p_bloblist->AddBlob(btp::Hadron_Decay);
  blob->AddToInParticles(inpart);
  SetPosition(blob);
  blob->SetTypeSpec("Sherpa");
  Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if (table==NULL) {
    msg_Error()<<METHOD<<" decay table not found, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  blob->AddData("dc",new Blob_Data<Decay_Channel*>(table->Select()));

  blob->AddData("p_onshell",new Blob_Data<Vec4D>(inpart->Momentum()));
  DEBUG_VAR(inpart->Momentum());
}

bool Hadron_Decay_Handler::
RejectExclusiveChannelsFromFragmentation(Blob* fblob)
{
  static std::string mname(METHOD);
  Blob * decayblob(NULL), * showerblob(NULL);
  if(fblob->Type()==btp::Fragmentation) {
    showerblob = fblob->UpstreamBlob();
    if(showerblob && showerblob->Type()==btp::Shower) {
      decayblob = showerblob->UpstreamBlob();
      if(decayblob && decayblob->Type()==btp::Hadron_Decay) {
        DEBUG_FUNC(decayblob->InParticle(0));
        DEBUG_VAR(*fblob);
        DEBUG_VAR(decayblob->InParticle(0)->Flav());
        Return_Value::IncCall(mname);
      }
    }
  }
  else if (fblob->Type()==btp::Hadron_Decay &&
	   (*fblob)["Partonic"]!=NULL) decayblob = fblob;
  if (decayblob) {
    bool anti=false;
    Decay_Map::iterator dt=
      p_decaymap->find(decayblob->InParticle(0)->Flav());
    if (dt==p_decaymap->end()) {
      anti=true;
      dt=p_decaymap->find(decayblob->InParticle(0)->Flav().Bar());
      if (dt==p_decaymap->end() || dt->second.size()!=1) {
	PRINT_INFO("Internal error.");
	throw Return_Value::Retry_Event;
      }
    }
    Flavour_Vector tmp, tmpno;
    for(int i=0;i<fblob->NOutP();i++) {
      Flavour flav(anti?fblob->OutParticle(i)->Flav().Bar():
		   fblob->OutParticle(i)->Flav());
      tmp.push_back(flav);
      if (!flav.IsPhoton()) tmpno.push_back(flav);
    }
    if (fblob!=decayblob) {
      for(int i=0;i<decayblob->NOutP();i++) {
	if (decayblob->OutParticle(i)->GetFlow(1)==0 &&
	    decayblob->OutParticle(i)->GetFlow(2)==0)
	  tmp.push_back(anti?decayblob->OutParticle(i)->Flav().Bar():
			decayblob->OutParticle(i)->Flav());
      }
    }
    std::sort(tmp.begin(), tmp.end(), Decay_Channel::FlavourSort);
    Flavour_Vector compflavs(tmp.size()+1), compflavsno(tmpno.size()+1);
    compflavs[0]=(anti?decayblob->InParticle(0)->Flav().Bar():
		  decayblob->InParticle(0)->Flav());
    for (size_t i=0; i<tmp.size(); ++i) compflavs[i+1]=tmp[i];
    if (tmpno.size()!=tmp.size()) {
      std::sort(tmpno.begin(), tmpno.end(), Decay_Channel::FlavourSort);
      compflavsno[0]=(anti?decayblob->InParticle(0)->Flav().Bar():
		      decayblob->InParticle(0)->Flav());
      for (size_t i=0; i<tmpno.size(); ++i) compflavsno[i+1]=tmpno[i];
    }
    Blob_Data_Base* data = (*fblob)["dc"];
    /*
      std::string decname(data->Get<Decay_Channel*>()->Name());
      bool outtag(decayblob->InParticle(0)->Flav()==Flavour(kf_B_s));
      if (outtag) {
      msg_Out()<<om::red
      <<METHOD<<" for "<<decayblob->InParticle(0)->Flav()<<" ->";
      for (size_t i=0; i<compflavs.size(); i++) 
      msg_Out()<<" "<<compflavs[i];
      msg_Out()<<"\n --> |"<<decname<<"|\n"
      <<(*fblob)<<"\n";
      }
    */
    if (dt->second[0]->GetDecayChannel(compflavs) ||
       (tmpno.size()!=tmp.size() && 
	dt->second[0]->GetDecayChannel(compflavsno))) {
      Return_Value::IncRetryPhase(mname);
      DEBUG_INFO("rejected. retrying decay.");
      if (fblob!=decayblob) p_bloblist->Delete(fblob);
      if (showerblob) p_bloblist->Delete(showerblob);
      decayblob->DeleteOutParticles();
      decayblob->InParticle(0)->SetStatus(part_status::active);
      Flavour infl(decayblob->InParticle(0)->Flav());
      decayblob->AddStatus(blob_status::internal_flag);
      Hadron_Decay_Table * hdt = 
	dynamic_cast<Hadron_Decay_Table*>(p_decaymap->FindDecay(infl));
      hdt->Select(decayblob);
      Spin_Density* sigma=
	m_spincorr ? new Spin_Density(decayblob->InParticle(0)) : NULL;
      Decay_Matrix* D=FillDecayTree(decayblob, sigma);
      delete sigma;
      delete D;
      /*
	if (outtag) {
	msg_Out()<<".\n"<<"Define new blob:"
	<<om::green<<(*decayblob)<<om::reset<<".\n"
	<<"=================================================\n";
	}
      */
      return true;
    }
    else {
      //if (outtag) msg_Out()<<om::green<<" --> accepted.\n"<<om::reset;
      DEBUG_INFO("not found as exclusive. accepted.");
      Vec4D vertex_position=decayblob->Position();
      if (showerblob) showerblob->SetPosition(vertex_position);
      if (fblob!=decayblob) fblob->SetPosition(vertex_position);
    }
  }
  return false;
}

void Hadron_Decay_Handler::SetPosition(ATOOLS::Blob* blob)
{
  Particle* inpart = blob->InParticle(0);
  if(inpart->Flav().Kfcode()==kf_K) {
    blob->SetPosition(inpart->XProd());
    return;
  }
  
  // boost lifetime into lab
  double gamma = 1./rpa->gen.Accu();
  if (inpart->Flav().HadMass()>rpa->gen.Accu()) {
    gamma = inpart->E()/inpart->Flav().HadMass();
  }
  else {
    double q2    = dabs(inpart->Momentum().Abs2());
    if (q2>rpa->gen.Accu()) gamma = inpart->E()/sqrt(q2);
  }
  double lifetime_boosted = gamma * inpart->Time();
  
  Vec3D      spatial = inpart->Distance( lifetime_boosted );
  Vec4D     position = Vec4D( lifetime_boosted*rpa->c(), spatial );
  blob->SetPosition( inpart->XProd() + position ); // in mm
}

