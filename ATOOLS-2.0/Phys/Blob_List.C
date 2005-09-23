#include "Blob_List.H"

#include "Blob.H"
#include "Particle.H"
#include "Message.H"

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
  m_destructor(NULL) {}

Blob_List::Blob_List(const bool destruct):
  m_destructor(destruct?this:NULL) {}

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
  if (deleted.find(blob)!=deleted.end()) return;
  deleted.insert(blob);
  for (int i=blob->NOutP()-1;i>=0;i=Min(blob->NOutP()-1,i-1)) {
    Blob *dblob=blob->OutParticle(i)->DecayBlob();
    if (dblob!=NULL) DeleteConnected(dblob,deleted);
  }
  for (int i=blob->NInP()-1;i>=0;i=Min(blob->NInP()-1,i-1)) {
    Blob *pblob=blob->InParticle(i)->ProductionBlob();
    if (pblob!=NULL) DeleteConnected(pblob,deleted);
  }
  delete blob;
  for (Blob_List::iterator bit=begin();bit!=end();++bit)
    if (*bit==blob) {
      erase(bit);
      break;
    }
}

size_t Blob_List::DeleteConnected(Blob *blob)
{
  if (blob==NULL) return 0;
  std::set<Blob*> deleted;
  for (Blob_List::iterator bit=begin();bit!=end();++bit) 
    if (*bit==blob) {
      DeleteConnected(blob,deleted);
      break;
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
	if (!TotalFourMomentum(part->ProductionBlob(),summed,inisum,finsum,mode))
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
    msg.Error()<<"Blob_List::FourMomentumConservation(): ("
	       <<this<<") Invalid momenta."<<std::endl;
    return false;
  }
  bool test=inisum==finsum;
  if (!test) {
    msg.Error()<<"Blob_List::FourMomentumConservation(): ("
	       <<this<<") Four Momentum is not conserved.\n"
	       <<"   p_{in} = "<<inisum<<" vs. p_{out} = "
	       <<finsum<<"."<<std::endl;
    if (msg.LevelIsDebugging()) {
      msg.Out()<<*this<<std::endl;
      for (Blob_List::const_iterator bit=begin();bit!=end();++bit) {
	Vec4D sum((*bit)->CheckMomentumConservation());
	if (sum!=Vec4D()) {
	  msg.Out()<<METHOD<<"(..): sum = "<<sum
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

void Blob_List::Clear() 
{
  while (!empty()) {
    delete back();
    pop_back();
  }
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
      if (flows.find(real)!=flows.end()) {
	msg.Error()<<"Blob_List::ColorConservation(): "
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
	msg.Error()<<"Blob_List::ColorConservation(): "
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
    msg.Error()<<"Blob_List::ColorConservation(): "
	       <<"Unconnected particles {\n";
    for (std::map<int,Particle*>::iterator uit=flows.begin();
	 uit!=flows.end();++uit) msg.Error()<<"   "<<*uit->second<<"\n";
    msg.Error()<<"}\n"<<*this<<std::endl;
    singlet=false;
  }
  return singlet;
}
