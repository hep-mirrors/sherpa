#include "ATOOLS/Phys/Blob_List.H"

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Smart_Pointer.H"

#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &s,const Blob_List &list) 
{
  s<<"Blob List with "<<list.size()<<" elements {"<<std::endl;
  {
    msg_Indent();
    for (Blob_List::const_iterator bit=list.begin(); bit!=list.end(); ++bit) {
      s<<*(*bit)<<std::endl;
    }
  }
  s<<"}"<<std::endl;
  return s;
}

Blob_List::Blob_List():
  m_destructor(NULL), m_extweight(1.) {}

Blob_List::Blob_List(const bool destruct):
  m_destructor(destruct?this:NULL), m_extweight(1.) {}

Blob *Blob_List::FindFirst(const btp::code code) const
{
  for (Blob_List::const_iterator bit=begin();bit!=end();++bit)
    if ((*bit)->Type()&code) return *bit;
  return NULL;
}

Blob *Blob_List::FindLast(const btp::code code) const
{
  for (Blob_List::const_reverse_iterator bit=rbegin();bit!=rend();++bit)
    if ((*bit)->Type()&code) return *bit;
  return NULL;
}

Blob_List Blob_List::Find(const btp::code code) const
{
  Blob_List hits;
  for (Blob_List::const_iterator bit=begin();bit!=end();++bit)
    if ((*bit)->Type()&code) hits.push_back(*bit);
  return hits;
}

void Blob_List::FindConnected(Blob *blob,Blob_List &connected,
			      std::set<const Blob*> &selected)
{
  if (selected.find(blob)!=selected.end()) return;
  selected.insert(blob);
  connected.push_back(blob);
  for (int i=blob->NOutP()-1;i>=0;i=Min(blob->NOutP()-1,i-1)) {
    Blob *dblob=blob->ConstOutParticle(i)->DecayBlob();
    if (dblob!=NULL) FindConnected(dblob,connected,selected);
  }
  for (int i=blob->NInP()-1;i>=0;i=Min(blob->NInP()-1,i-1)) {
    Blob *pblob=blob->ConstInParticle(i)->ProductionBlob();
    if (pblob!=NULL) FindConnected(pblob,connected,selected);
  }
}

Blob_List Blob_List::FindConnected(const Blob *blob)
{
  Blob_List connected;
  if (blob==NULL) return connected;
  std::set<const Blob*> selected;
  FindConnected((Blob*)blob,connected,selected);
  return connected;
}

Blob_List Blob_List::FindConnected(const Particle *particle)
{
  if (particle==NULL) return Blob_List();
  Blob *owner=particle->DecayBlob();
  if (owner==NULL) owner=particle->ProductionBlob();
  if (owner==NULL) return Blob_List();
  return FindConnected(owner);
}

bool Blob_List::Delete(Blob *blob) 
{
  if (blob==NULL) return false;
  for (Blob_List::iterator bit=begin();bit!=end();++bit) 
    if (*bit==blob) {
      erase(bit);
      blob->RemoveOwnedParticles();
      delete blob;
      return true;
    }
  return false;
}

void Blob_List::DeleteConnected(Blob *blob,std::set<Blob*> &deleted)
{
  if (blob==NULL || deleted.find(blob)!=deleted.end()) return;
  deleted.insert(blob);
  Particle_Vector parts(blob->GetInParticles());
  for (Particle_Vector::iterator pit(parts.begin());pit!=parts.end();++pit) 
    DeleteConnected((*pit)->ProductionBlob(),deleted);
  parts=blob->GetOutParticles();
  for (Particle_Vector::iterator pit(parts.begin());pit!=parts.end();++pit) 
    DeleteConnected((*pit)->DecayBlob(),deleted);
}

size_t Blob_List::DeleteConnected(Blob *blob)
{
  std::set<Blob*> deleted;
  DeleteConnected(blob,deleted);
  for (Blob_List::iterator bit(begin());bit!=end();++bit) {
    std::set<Blob*>::const_iterator rit(deleted.find(*bit));
    if (rit!=deleted.end()) {
      delete *bit;
      --(bit=erase(bit));
    }
  }
  return deleted.size();
}

size_t Blob_List::DeleteConnected(Particle *particle)
{
  if (particle==NULL) return 0;
  Blob *owner=particle->DecayBlob();
  if (owner==NULL) owner=particle->ProductionBlob();
  if (owner==NULL) return 0;
  return DeleteConnected(owner);
}

bool Blob_List::TotalFourMomentum(Blob *blob,std::set<Blob*> &summed,
				  Vec4D &inisum,Vec4D &finsum,
				  const int mode) const
{
  if (summed.find(blob)!=summed.end()) return true;
  summed.insert(blob);
  bool success=true;
  if (mode<=0)
    for (int i=0;i<blob->NInP();++i) {
      const ATOOLS::Particle *part=blob->ConstInParticle(i);
      double abs2=part->Momentum().Abs2();
      if (abs2>0 && abs2<0) return false;
      if (part->ProductionBlob()==NULL) inisum+=part->Momentum(); 
      else 
	if (!TotalFourMomentum(part->ProductionBlob(),
			       summed,inisum,finsum,mode))
	  success=false;
    }
  if (mode>=0)
    for (int i=0;i<blob->NOutP();++i) {
      const ATOOLS::Particle *part=blob->ConstOutParticle(i);
      double abs2=part->Momentum().Abs2();
      if (abs2>0 && abs2<0) return false;
      if (part->DecayBlob()==NULL) finsum+=part->Momentum(); 
      else 
	if (!TotalFourMomentum(part->DecayBlob(),summed,inisum,finsum,mode))
	  success=false;
    }
  return success;
}

Vec4D Blob_List::TotalFourMomentum() const
{
  if (empty()) return Vec4D();
  Vec4D inisum,finsum;
  std::set<ATOOLS::Blob*> summed;
  if (!TotalFourMomentum(*begin(),summed,inisum,finsum,0)) 
    return Vec4D(sqrt(-1.0),Vec3D());
  return finsum-inisum;
}

Vec4D Blob_List::IncomingFourMomentum() const
{
  if (empty()) return Vec4D();
  Vec4D inisum,finsum;
  std::set<ATOOLS::Blob*> summed;
  if (!TotalFourMomentum(*begin(),summed,inisum,finsum,-1))
    return Vec4D(sqrt(-1.0),Vec3D());
  return inisum;
}

Vec4D Blob_List::OutgoingFourMomentum() const
{
  if (empty()) return Vec4D();
  Vec4D inisum,finsum;
  std::set<ATOOLS::Blob*> summed;
  if (!TotalFourMomentum(*begin(),summed,inisum,finsum,1))
    return Vec4D(sqrt(-1.0),Vec3D());
  return finsum;
}

bool Blob_List::FourMomentumConservation() const
{
  if (empty()) return true;
  Vec4D inisum,finsum;
  std::set<ATOOLS::Blob*> summed;
  if (!TotalFourMomentum(*begin(),summed,inisum,finsum,0)) {
    msg_Error()<<METHOD<<"(): ("<<this<<") Invalid momenta."<<std::endl;
    return false;
  }
  static double accu(sqrt(rpa->gen.Accu()));
  bool test=IsEqual(inisum,finsum,accu);
  if (!test) {
    msg_Error()<<METHOD<<"(): ("<<this<<") Four Momentum is not conserved.\n"
               <<"   p_{in} = "<<inisum<<" vs. p_{out} = "<<finsum<<",\n"
               <<"   diff = "<<finsum-inisum<<"."<<std::endl;
    if (msg_LevelIsDebugging()) {
      msg_Out()<<*this<<std::endl;
      for (Blob_List::const_iterator bit=begin();bit!=end();++bit) {
	Vec4D sum((*bit)->CheckMomentumConservation());
	if (sum!=Vec4D()) {
	  msg_Out()<<METHOD<<"(..): sum = "<<sum
		   <<" in\n"<<**bit<<std::endl;
	}
      }
    }
  }
  return test;
}

Particle_List Blob_List::ExtractParticles(const int status,
					  const int mode) const
{
  Particle_List particles;
  for (Blob_List::const_iterator bit=begin();bit!=end();++bit) {
    if (mode>=0)
      for (int i=0;i<(*bit)->NOutP();++i) {
	ATOOLS::Particle *part=(*bit)->OutParticle(i);
	if (part->Status()==status) particles.push_back(part);
      }
    if (mode<=0)
      for (int i=0;i<(*bit)->NInP();++i) {
	ATOOLS::Particle *part=(*bit)->InParticle(i);
	if (part->Status()==status) particles.push_back(part);
      }
  }
  return particles;
}

Particle_List Blob_List::ExtractLooseParticles(const int mode) const
{
  Particle_List particles;
  for (Blob_List::const_iterator bit=begin();bit!=end();++bit) {
    if (mode>=0)
      for (int i=0;i<(*bit)->NOutP();++i) {
	ATOOLS::Particle *part=(*bit)->OutParticle(i);
	if (part->DecayBlob()==NULL) particles.push_back(part);
      }
    if (mode<=0)
      for (int i=0;i<(*bit)->NInP();++i) {
	ATOOLS::Particle *part=(*bit)->InParticle(i);
	if (part->ProductionBlob()==NULL) particles.push_back(part);
      }
  }
  return particles;
}

void Blob_List::Clear(Blob * blob) 
{
  if (blob==NULL) {
    while (!empty()) {
      if (back()) delete back();
      pop_back();
    }
    return;
  }
  for (int i(0);i<blob->NInP();++i) 
    if (blob->InParticle(i)->ProductionBlob()!=NULL)
      blob->InParticle(i)->ProductionBlob()->
	RemoveOutParticle(blob->InParticle(i));
  for (int i(0);i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->DecayBlob()!=NULL)
      blob->OutParticle(i)->DecayBlob()->
	RemoveInParticle(blob->OutParticle(i));
    blob->OutParticle(i)->SetStatus(part_status::active);
  }
  for (const_iterator bit(begin());bit!=end();++bit) 
    if (*bit!=blob) delete *bit;
  resize(1);
  back()=blob;
}

bool Blob_List::ColorConservation() const
{
  bool singlet=true;
  Particle_List outgoing=ExtractLooseParticles();
  std::map<int,Particle*> flows;
  for (Particle_List::const_iterator pit=outgoing.begin();
       pit!=outgoing.end();++pit) {
    int real=(*pit)->GetFlow()->Code(1);
    int anti=-(*pit)->GetFlow()->Code(2);
    if (real!=0) {
      if (anti!=0 && real==-anti) {
	msg_Error()<<"Blob_List::ColorConservation(): "
		   <<"Color singlet gluon "<<**pit<<std::endl;
	return false;
      }
      if (flows.find(real)!=flows.end()) {
	msg_Error()<<"Blob_List::ColorConservation(): "
			   <<"Doubled color index '"<<real<<"' {\n   "
			   <<**pit<<"\n   "<<*flows[real]<<"\n}"<<std::endl;
	singlet=false;
      }
      std::map<int,Particle*>::iterator dit=flows.find(-real);
      if (dit!=flows.end()) flows.erase(dit);
      else flows[real]=*pit;
    }
    if (anti!=0) {
      if (flows.find(anti)!=flows.end()) {
	msg_Error()<<"Blob_List::ColorConservation(): "
			   <<"Doubled color index '"<<anti<<"' {\n   "
			   <<**pit<<"\n   "<<*flows[anti]<<"\n}"<<std::endl;
	singlet=false;
      }
      std::map<int,Particle*>::iterator dit=flows.find(-anti);
      if (dit!=flows.end()) flows.erase(dit);
      else flows[anti]=*pit;
    }
  }
  if (!flows.empty()) {
    msg_Error()<<"Blob_List::ColorConservation(): "
	       <<"Unconnected particles {\n";
    for (std::map<int,Particle*>::iterator uit=flows.begin();
	 uit!=flows.end();++uit) msg_Error()<<"   "<<*uit->second<<"\n";
    msg_Error()<<"}\n"<<*this<<std::endl;
    singlet=false;
  }
  return singlet;
}

Blob *Blob_List::AddBlob(const btp::code &type)
{
  Blob *blob(new Blob());
  blob->SetType(type);
  blob->SetId();
  blob->SetStatus(blob_status::inactive);
  push_back(blob);
  return blob;
}

bool Blob_List::MergeSubsequentType(btp::code mtype,btp::code dtype,
				    long int & NBlob, long int & NPart) {
  bool merger(false);
  Blob_List::iterator mother(begin()),daughter;
  while (mother!=end()) {
    if ((*mother)->Type()==mtype) {
      for (int i=0;i<(*mother)->NOutP();i++) {
	Particle * part((*mother)->OutParticle(i));
	Blob * blob(part->DecayBlob());
	if (blob && blob->Type()==dtype) {
	  merger=true;
	  while (blob->NOutP()>0) {
	    (*mother)->AddToOutParticles(blob->RemoveOutParticle(blob->NOutP()-1,true));
	  }
	  daughter=begin();
	  while (daughter!=end()) {
	    if ((*daughter)==blob) {
	      NBlob--;
	      delete blob; 
	      daughter = erase(daughter);
	      break;
	    }
	    else daughter++;
	  }
	  NPart--;
	  (*mother)->DeleteOutParticle(part);
	}
      }
    }
    mother++;
  }
  return merger;
}

void Blob_List::MergeSubsequentTypeRecursively(btp::code mtype,btp::code dtype,
					       long int & NBlob, long int & NPart) {
  while (MergeSubsequentType(mtype,dtype,NBlob,NPart)) {}
}

bool Blob_List::Shorten(kf_code kfc,btp::code out,btp::code in)
{
  bool shorten(false);
  Blob_List::iterator mother(begin()),daughter;
  while (mother!=end()) {
    if (out==btp::Unspecified || (*mother)->Type()==out) {
      for (int i=0;i<(*mother)->NOutP();i++) {
	Particle * part((*mother)->OutParticle(i));
	if (part->Flav().Kfcode()!=kfc) continue;
	Blob * blob(part->DecayBlob());
	if (blob && (in==btp::Unspecified || blob->Type()==in)) {
	  shorten=true;
	  while (blob->NOutP()>0) {
	    (*mother)->AddToOutParticles(blob->RemoveOutParticle(blob->NOutP()-1,true));
	  }
	  daughter=begin();
	  while (daughter!=end()) {
	    if ((*daughter)==blob) {
	      delete blob; 
	      daughter = erase(daughter);
	      break;
	    }
	    else daughter++;
	  }
	  (*mother)->DeleteOutParticle(part);
	}
      }
    }
    mother++;
  }  
  return shorten;
}

void Blob_List::ShortenRecursively(kf_code kfc,btp::code out,btp::code in) {
  while (Shorten(kfc,out,in)) {}
}

double Blob_List::Weight() const
{
  double weight(1.0);
  bool noweight(true);
  for (const_iterator it(begin());it!=end();++it) {
    Blob_Data_Base *bd((**it)["Weight"]);
    if (bd) {
      weight*=bd->Get<double>();
      noweight=false;
    }
  }
  return noweight?m_extweight:weight;
}

Blob_List Blob_List::Copy() const
{
  Blob_List copy;
  copy.resize(size());
  std::map<Particle*,Particle*> pmap;
  for (size_t i(0);i<size();++i) {
    Blob *cb((*this)[i]), *nb(copy[i] = new Blob(cb,false));
    for (int j(0);j<cb->NInP();++j) {
      Particle *cp(cb->InParticle(j)), *np(NULL);
      std::map<Particle*,Particle*>::iterator pit(pmap.find(cp));
      if (pit!=pmap.end()) np=pit->second;
      else pmap[cp]=np = new Particle(*cp);
      if (np->DecayBlob()) THROW(fatal_error,"Internal error");
      nb->AddToInParticles(np);
    }
    for (int j(0);j<cb->NOutP();++j) {
      Particle *cp(cb->OutParticle(j)), *np(NULL);
      std::map<Particle*,Particle*>::iterator pit(pmap.find(cp));
      if (pit!=pmap.end()) np=pit->second;
      else pmap[cp]=np = new Particle(*cp);
      if (np->ProductionBlob()) THROW(fatal_error,"Internal error");
      nb->AddToOutParticles(np);
    }
  }

  // adjust p_originalpart pointers after copying
  for (size_t i=0; i<copy.size(); ++i) {
    Blob* nb=copy[i];
    for (int j=0;j<nb->NInP();++j) {
      std::map<Particle*,Particle*>::iterator pit(pmap.find(nb->InParticle(j)->OriginalPart()));
      if (pit!=pmap.end()) nb->InParticle(j)->SetOriginalPart(pit->second);
      else nb->InParticle(j)->SetOriginalPart(nb->InParticle(j));
    }
    for (int j=0;j<nb->NOutP();++j) {
      std::map<Particle*,Particle*>::iterator pit(pmap.find(nb->OutParticle(j)->OriginalPart()));
      if (pit!=pmap.end()) nb->OutParticle(j)->SetOriginalPart(pit->second);
      else nb->OutParticle(j)->SetOriginalPart(nb->OutParticle(j));
    }
  }

  // adjust particle pointers in amplitude tensor
  Blob* signal=copy.FindFirst(btp::Signal_Process);
  if (signal) {
    Blob_Data_Base* data = (*signal)["ATensor"];
    if (data) {
      SP(METOOLS::Amplitude2_Tensor) origamps = data->Get<SP(METOOLS::Amplitude2_Tensor)>();
      SP(METOOLS::Amplitude2_Tensor) newamps(new METOOLS::Amplitude2_Tensor(*origamps));
      newamps->UpdateParticlePointers(pmap);
      data->Set(newamps);
    }
  }
  return copy;
}
