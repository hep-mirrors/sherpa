#include"Hard_Decay_Handler.H"

#include "Full_Decay_Table.H"
#include "Data_Read.H"
#include "Message.H"
#include "MyStrStream.H"

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace SHERPA;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Hard_Decay_Handler::Hard_Decay_Handler(std::string _path,std::string _file,std::string _pfile,
				       MODEL::Model_Base * _model) :
   m_on(0), m_path(_path), m_file(_file), p_mehandler(NULL), m_meflag(0), p_amegic(NULL)
{
  ReadInDecays();

  if (m_decaytables.size()>0) {
    m_on = 1;
    EvaluateWidths(_pfile,_model);
    SetWidths();
  }
}

Hard_Decay_Handler::~Hard_Decay_Handler() { }

void Hard_Decay_Handler::ReadInDecays()
{
  ifstream from((m_path+m_file).c_str());
  if (!from) {
    msg.Out()<<"Warning :  in Hard_Decay_Handler::ReadInDecays : "<<endl
	     <<"   File : "<<(m_path+m_file)<<" not found ! "<<endl;
    return;
  }

  char   buffer[100];
  string buf,number;
  int    kfc;
  size_t pos;
  Decay_Table * dt = NULL;
  Flavour     flav;
  FlavourSet  decflavs;
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '%' && strlen(buffer)>0) {
      buf    = string(buffer);
      // Init decay table for another particle
      pos    = buf.find(string("Decays :")); 
      if (pos!=std::string::npos && pos<=buf.length()) {
	buf  = buf.substr(pos+8);
	while(buf.length()>0) {
	  if (buf[0]==' ') buf = buf.substr(1);
	  else {
	    pos = buf.find(string(" "));
	    if (pos!=std::string::npos) buf = buf.substr(0,pos);
	    MyStrStream sstream;
	    sstream<<buf;
	    sstream>>kfc;
	    break;
	  }
	}
	flav = Flavour(kf::code(int(abs(double(kfc)))));
	dt   = new Decay_Table(Flavour(kf::code(int(abs(double(kfc))))));
	m_decaytables.insert(dt);
      }

      pos     = buf.find(string("forced channel :"));  
      if (pos!=std::string::npos && pos<=buf.length()) {
	decflavs.clear();
	buf  = buf.substr(pos+16);
	buf  = buf+string(" ");
	while(buf.length()>0) {
	  if (buf[0]==' ') buf = buf.substr(1);
	  else {
	    pos = buf.find(string(" "));
	    MyStrStream sstream;
	    sstream<<buf.substr(0,pos);
	    sstream>>kfc;
	    flav = Flavour(kf::code(int(abs(double(kfc)))));
	    if (kfc<0) flav = flav.Bar();
	    decflavs.insert(flav);
	    buf = buf.substr(pos);
	  }
	}
	if (dt) dt->SetSelectedChannel(decflavs);
      }

      // Check, if width of particle.dat is to be overwritten by total width as calculated
      pos     = buf.find(string("overwrite"));  
      if (pos!=std::string::npos && pos<=buf.length()) dt->SetOverwrite(); 

      // Check, if Breit Wigner smearing is to be applied in event generation
      pos     = buf.find(string("Breit-Wigner on"));  
      if (pos!=std::string::npos && pos<=buf.length()) dt->SetSmearing(); 
    }
  }
}

void Hard_Decay_Handler::EvaluateWidths(std::string _pfile,MODEL::Model_Base * _model)
{
  Flavour    flav;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if ((*dit)->Overwrite()) {
      if (_model->FillDecay((*dit))) { 
	if (msg.LevelIsTracking()) { (*dit)->Output(); }
      }
      else {
	(*dit)->Flav().SetWidth(-1.);
	if (!p_mehandler) p_mehandler = new Matrix_Element_Handler(m_path,_pfile,_model,NULL);
	flav = (*dit)->Flav();
	if (!p_mehandler->AddToDecays(flav)) {
	  msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		     <<"   Could not add "<<flav
		     <<" to list of decays treated by ME_Handler. Abort run."<<endl;
	  abort();
	}
      }
    }
    else if ((*dit)->FixedDecay()) {
      if (!p_mehandler) p_mehandler = new Matrix_Element_Handler(m_path,_pfile,_model,NULL);
      if (!p_mehandler->AddToDecays((*dit)->GetOneDecayChannel())) {
	msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		   <<"   Could not add "<<flav
		   <<" to list of decays treated by ME_Handler. Abort run."<<endl;
	abort();
      }
    }
    else {
      if (!p_mehandler) p_mehandler = new Matrix_Element_Handler(m_path,_pfile,_model,NULL);
      flav = (*dit)->Flav();
      if (!p_mehandler->AddToDecays(flav)) {
	msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		   <<"   Could not add "<<flav
		   <<" to list of decays treated by ME_Handler. Abort run."<<endl;
	abort();
      }
    }
  }
  if (p_mehandler->InitializeDecayTables()) p_mehandler->CalculateWidths();
}


void Hard_Decay_Handler::SetWidths()
{
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    p_mehandler->FillDecayTable((*dit),(*dit)->Overwrite());
    //if (!flag && (*dit)->Overwrite()) p_mehandler->FillDecayTable((*dit),true);
    //if (flag && !(*dit)->Overwrite()) p_mehandler->FillDecayTable((*dit),false);
  }
} 

double Hard_Decay_Handler::DefineSecondaryDecays(ATOOLS::Blob * _blob,bool _add) 
{
  Particle      * particle;
  Decay_Channel * dc = NULL;
  double Mmin = 0.;
  double rest = (_add)? 0. : _blob->InParticle(0)->Flav().Mass();
  DecayingParticleTable dptable;
  for (short int i=0;i<_blob->NOutP();i++) {
    particle       = _blob->OutParticle(i);
    if (_add) rest += particle->Momentum()[0];
    if (particle->Flav().IsStable()) Mmin+=particle->Flav().Mass();
                                else dptable.insert(make_pair(particle,dc));
  }
  rest       -= Mmin; 
  for (DPTIt dptit=dptable.begin();dptit!=dptable.end();++dptit) {
    dptit->second = SpecifyHardDecay(dptit->first,rest);
    m_table.insert(make_pair(dptit->first,dptit->second));
  }
  bool smearit;
  for (DPTIt dptit=dptable.begin();dptit!=dptable.end();++dptit) {
    smearit = false;
    for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
      if ((*dit)->Flav()==dptit->first->Flav() ||
	  (*dit)->Flav().Bar()==dptit->first->Flav()) {
	if ((*dit)->Smearing()) smearit = true;
      }
    }
    if (smearit) dptit->first->SetFinalMass(dptit->first->FinalMass(),rest);
            else dptit->first->SetFinalMass(dptit->first->Flav().Mass());
    Mmin += dptit->first->FinalMass();
    rest -= dptit->first->FinalMass()-dptit->second->MinimalMass();
  }  
  if (_add) {
    _blob->BoostInCMS();
    p_tools->ShuffleMomenta(_blob->GetOutParticles());
    _blob->BoostInLab();
  }
  return Mmin;
} 

bool Hard_Decay_Handler::PerformDecay(ATOOLS::Blob * _blob) {
  Decay_Channel * dc=0;
  Particle      * part = _blob->InParticle(0);
  for (DPTIt dptit=m_table.begin();dptit!=m_table.end();++dptit) {
    if (dptit->first->Flav()==part->Flav() &&
	dptit->first->FinalMass()==part->FinalMass()) {
      dc = dptit->second;
      break;
    }
  }
  p_mehandler->GenerateOneEvent(dc,part->FinalMass());
  Poincare lab(_blob->InParticle(0)->Momentum());

  bool shuffle = false;
  Vec4D mom;
  for (unsigned int i=0;i<p_mehandler->NDecOut();i++) {
    mom = p_mehandler->DecMomenta()[i+1];
    lab.BoostBack(mom);
    _blob->OutParticle(i)->SetMomentum(mom);
    if (!_blob->OutParticle(i)->Flav().IsStable()) shuffle = true;
  }
  if (shuffle) p_tools->ShuffleMomenta(_blob->GetOutParticles());
  return true;
}

Decay_Channel * Hard_Decay_Handler::SpecifyHardDecay(ATOOLS::Particle * _part,double & _mmax) 
{
  bool barflag  = false;
  bool unstable = false;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if ((*dit)->Flav()==_part->Flav() ||
	(*dit)->Flav().Bar()==_part->Flav()) {
      if ((*dit)->Flav().Bar()==_part->Flav() &&
	  (*dit)->Flav()!=_part->Flav()) barflag = true;
      (*dit)->Select();
      Decay_Channel * dc = (*dit)->GetOneDecayChannel();
      Blob * blob         = new Blob();
      blob->AddToInParticles(_part);
      blob->SetType(btp::Hard_Decay);
      blob->SetTypeSpec(dc->ProcessName());
      blob->SetBeam(-1);
      blob->SetStatus(1);
      _part->SetDecayBlob(blob);
      Particle * particle;
      Flavour flav;
      for (int i=0;i<dc->NumberOfDecayProducts();i++) {
	particle = new Particle();
	flav     = dc->GetDecayProduct(i);
	if (barflag) flav = flav.Bar();
	if (!flav.IsStable()) unstable = true;
	particle->SetFlav(flav);
	particle->SetInfo('H');
	particle->SetProductionBlob(blob);
	particle->SetNumber((long int)particle);
	blob->AddToOutParticles(particle);
      }
      double decmass;
      if (unstable) decmass = DefineSecondaryDecays(blob,0);
               else decmass = dc->MinimalMass();

      _part->SetFinalMass(decmass);
      _mmax          -= decmass;
      return dc;
    }
  }

  msg.Error()<<"Error in Hard_Decay_Handler::SpecifyHardDecay("<<_part->Flav()<<")"<<endl
	     <<"   No Decay_Table found in list. Abort."<<endl;
  abort();
}

std::string Hard_Decay_Handler::Name() 
{
  return std::string("");
}

void Hard_Decay_Handler::ResetTables()
{
  m_table.clear();
} 

