#include "Beam_Remnant_Handler.H"

#include "Hadron_Remnant.H"
#include "Electron_Remnant.H"
#include "Photon_Remnant.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace ATOOLS;

#ifdef DEBUG__Beam_Remnant_Handler
std::set<Blob*> checked_blobs;

void SumMomenta(Blob * bl, Vec4D & inisum, Vec4D & finsum,bool iterate=true) 
{
  // check if caught in loop
  std::set<Blob*>::const_iterator bit=checked_blobs.find(bl);
  if (bit!=checked_blobs.end()) return;
  checked_blobs.insert(bl);
  for (int i=0;i<bl->NInP();++i) {
    Particle * p =bl->InParticle(i);
    if (p->ProductionBlob()==NULL) inisum+=p->Momentum();
    else if (iterate) SumMomenta(p->ProductionBlob(),inisum,finsum);
    else inisum+=p->Momentum();
  }
  for (int i=0;i<bl->NOutP();++i) {
    Particle * p =bl->OutParticle(i);
    if (p->DecayBlob()==NULL) finsum+=p->Momentum();
    else if (iterate) SumMomenta(p->DecayBlob(),inisum,finsum);
    else finsum+=p->Momentum();
  }
}

bool SumMomenta(Blob * bl)
{
  Vec4D inisum,finsum;
  checked_blobs.clear();
  SumMomenta(bl,inisum,finsum);
  std::cout<<" inisum="<<inisum<<std::endl;
  std::cout<<" finsum="<<finsum<<std::endl;
  return (inisum==finsum);
}
#endif

Beam_Remnant_Handler::Beam_Remnant_Handler(std::string _m_path,std::string _m_file,
					   PDF::ISR_Handler * _p_isr,BEAM::Beam_Spectra_Handler * _p_beam,
					   double _m_scale) :
  p_isr(_p_isr), 
  p_beam(_p_beam),
  m_path(_m_path), 
  m_file(_m_file),
  m_fill(true)
{
  for (size_t i=0;i<2;++i) {
    if (p_isr->Flav(i).IsHadron()) p_beampart[i] = new Hadron_Remnant(p_isr,i,_m_scale);
    else if (p_isr->Flav(i).IsLepton()) p_beampart[i] = new Electron_Remnant(p_isr,i,_m_scale);
    else if (p_isr->Flav(i).IsPhoton()) p_beampart[i] = new Photon_Remnant(i);
    else {
      ATOOLS::msg.Error()<<"Beam_Remnant_Handler::Beam_Remnant_Handler(..): "
			 <<"Cannot determine type of beam "<<i<<"! Abort."<<std::endl;
      exit(129);
    }
  }
  for (size_t i=0;i<2;++i) p_beampart[i]->SetPartner(p_beampart[1-i]);
  p_kperp = new Primordial_KPerp(_m_path,_m_file);
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  for (size_t i=0;i<2;++i) delete p_beampart[i];
  delete p_kperp;
}


bool Beam_Remnant_Handler::FillBunchBlobs(Blob_List * _bloblist,Particle_List * _particlelist)
{
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool flag=false;
  int  number;
  Particle * p;
  for (int i=0;i<2;i++) {
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      if ((*biter)->Status()==1 && (*biter)->Beam()==i && 
	  ((*biter)->Type()==btp::Beam || (*biter)->Type()==btp::IS_Shower)) {
	(*biter)->SetStatus(2);
	blob = new ATOOLS::Blob();
	blob->SetId(_bloblist->size());
	blob->SetType(btp::Bunch);
	blob->SetBeam(i);
	blob->AddToOutParticles((*biter)->InParticle(0));
	if ((*biter)->InParticle(0)->Flav()==p_beam->GetBeam(i)->Beam() &&
	    IsEqual((*biter)->InParticle(0)->E(),p_beam->GetBeam(i)->InMomentum()[0])) {
	  p = new Particle((*biter)->InParticle(0));
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);
	  blob->AddToInParticles(p);
	}
	else {
	  p = new Particle(-1,p_beam->GetBeam(i)->Beam(),p_beam->GetBeam(i)->InMomentum());
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);	  
	  p->SetStatus(2);
	  blob->AddToInParticles(p);
	  Particle * p = new Particle(-1,p_beam->GetBeam(i)->Remnant(),
				      p_beam->GetBeam(i)->InMomentum()+(-1.)*(*biter)->InParticle(0)->Momentum());
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);
	  blob->AddToOutParticles(p);
	}
	_bloblist->insert(_bloblist->begin(),blob);
	flag=true;
      }
    }
  }
#ifdef DEBUG__Beam_Remnant_Handler
  SumMomenta(_bloblist->front());
#endif
  return flag;
}

bool Beam_Remnant_Handler::FillBeamBlobs(Blob_List * _bloblist,Particle_List * _particlelist)
{ 
  if (!m_fill) return false;
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool okay=false;
  int number;
  for (int i=0;i<2;i++) {
    p_beampart[i]->Clear();
    bool flag=true,treat=false;
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      if ((*biter)->Beam()==i && (*biter)->Status()==1 && 
	  (*biter)->Type()==btp::IS_Shower) { 
	blob=NULL;
	if (p_beampart[i]->Type()==Remnant_Base::Hadron) {
	  if (flag) {
	    blob = new ATOOLS::Blob();
	    blob->SetId(_bloblist->size());
	    blob->SetType(btp::Beam);
	    blob->SetBeam(i);
	    blob->SetStatus(1);
	    Particle *p = new Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	    p->SetStatus(2);
	    blob->AddToInParticles(p);
	    _bloblist->insert(_bloblist->begin(),blob);
	    p_beamblob[i]=blob;
	    flag=false;
	    okay=true;
	  }
	  p_beampart[i]->ExtractParton((*biter)->InParticle(0));
	  (*biter)->SetStatus(2);
	  treat=true;
	}
	else if (!IsEqual((*biter)->InParticle(0)->E(),p_beam->GetBeam(i)->OutMomentum()[0])) {
	  (*biter)->SetStatus(2);
	  blob = new ATOOLS::Blob();
	  blob->SetId(_bloblist->size());
	  blob->SetType(btp::Beam);
	  blob->SetBeam(i);
	  blob->SetStatus(1);
	  p_beampart[i]->ExtractParton((*biter)->InParticle(0));
	  Particle *p = new Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);
	  p->SetStatus(2);
	  blob->AddToInParticles(p);
	  _bloblist->insert(_bloblist->begin(),blob);
	  p_beamblob[i]=blob;
	  flag=false;
	  okay=true;
	}
      }
    }
    if (blob!=NULL) p_beampart[i]->FillBlob(blob,_particlelist);
  }
  if ((p_beampart[0]->Type()==Remnant_Base::Hadron)||
      (p_beampart[1]->Type()==Remnant_Base::Hadron)) {
    p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
    for (size_t i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
  }
  for (size_t i=0;i<2;++i) p_beampart[i]->AdjustKinematics();
#ifdef DEBUG__Beam_Remnant_Handler
  SumMomenta(_bloblist->front());
#endif
  return okay;
}

