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
    SetWidths(0);
  }
}

Hard_Decay_Handler::~Hard_Decay_Handler() { }

void Hard_Decay_Handler::ReadInDecays()
{
  ifstream from((m_path+m_file).c_str());
  if (!from) {
    msg.Error()<<"Error in Hard_Decay_Handler::ReadInDecays : "<<endl
	       <<"   File : "<<(m_path+m_file).c_str()<<" not found ! Abort program execution."<<endl;
    abort();
  }

  char buffer[100];
  string buf,number;
  int pos,kfc;
  Decay_Table * dt = NULL;
  Flavour flav;
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '%' && strlen(buffer)>0) {
      buf    = string(buffer);
      pos    = buf.find(string("Decays :")); 
      if (pos>-1 && pos<(int)buf.length()) {
	buf  = buf.substr(pos+8);
	while(buf.length()>0) {
	  if (buf[0]==' ') buf = buf.substr(1);
	  else {
	    pos = buf.find(string(" "));
	    if (pos>0) buf = buf.substr(0,pos);
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
      pos     = buf.find(string("overwrite"));  
      if (pos>-1 && pos<(int)buf.length() && (dt)) dt->SetOverwrite(); 
    }
  }
}

void Hard_Decay_Handler::EvaluateWidths(std::string _pfile,MODEL::Model_Base * _model)
{
  Flavour flav;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if ((*dit)->Overwrite()) {
      if (_model->FillDecay((*dit))) { if (rpa.gen.Tracking()) { (*dit)->Output(); }}
      else {
	(*dit)->Flav().SetWidth(-1.);
	if (!p_mehandler) {
	  p_mehandler = new Matrix_Element_Handler(m_path,_pfile,_model,NULL);
	}
	flav = (*dit)->Flav();
	if (!p_mehandler->AddToDecays(flav)) {
	  msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		     <<"   Could not add "<<flav
		     <<" to list of decays treated by ME_Handler. Abort run."<<endl;
	  abort();
	}
      }
    }
  }
  if (p_mehandler->InitializeDecayTables()) p_mehandler->CalculateWidths();
}


void Hard_Decay_Handler::SetWidths(bool flag)
{
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if (!flag) {
      if ((*dit)->Overwrite()) {
	if ((*dit)->Flav().Width()<0.) p_mehandler->FillDecayTable((*dit),true);
      }
    }
    if (flag && !(*dit)->Overwrite()) p_mehandler->FillDecayTable((*dit),false);
    //if (rpa.gen.Tracking()) (*dit)->Output();
  }
} 

bool Hard_Decay_Handler::InitializeAllHardDecays(std::string _pfile,MODEL::Model_Base * _model) 
{
  bool newones = 0;
  Flavour flav;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    flav = (*dit)->Flav();
    if (p_mehandler->AddToDecays(flav)) newones = 1;
  }
  if (newones) {
    p_mehandler->InitializeDecayTables();
    p_mehandler->CalculateWidths();
  }
  SetWidths(1);

  return 1;
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
  for (DPTIt dptit=dptable.begin();dptit!=dptable.end();++dptit) {
    dptit->first->SetFinalMass(dptit->first->FinalMass(),rest);
    Mmin += dptit->first->FinalMass();
    rest -= dptit->first->FinalMass()-dptit->second->MinimalMass();
  }  
  if (_add) p_tools->ShuffleMomenta(_blob->GetOutParticles());
  return Mmin;
} 

bool Hard_Decay_Handler::PerformDecay(ATOOLS::Blob * _blob) {
  Decay_Channel * dc;
  Particle      * part = _blob->InParticle(0);
  for (DPTIt dptit=m_table.begin();dptit!=m_table.end();++dptit) {
    if (dptit->first->Flav()==part->Flav() &&
	dptit->first->FinalMass()==part->FinalMass()) {
      dc = dptit->second;
      break;
    }
  }
  //if (rpa.gen.Tracking()) dc->Output();
  p_mehandler->GenerateOneEvent(dc,part->FinalMass());
  bool shuffle = false;
  for (unsigned int i=0;i<p_mehandler->NDecOut();i++) {
    _blob->OutParticle(i)->SetMomentum(p_mehandler->DecMomenta()[i+1]);
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
      blob->SetType(string("Hard decay : ")+dc->ProcessName());
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
	particle->SetProductionBlob(blob);
	particle->SetNumber(int(particle));
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

