#include "Blob.H"
#include "Particle.H"
#include "Poincare.H"
#include "Message.H"
#include <iomanip>

using namespace ATOOLS;

std::ostream& ATOOLS::operator<<(std::ostream& ostr, const btp::code btpc) {
  switch (btpc) {
  case btp::Unspecified:        return ostr<<"Unspecified       ";
  case btp::Signal_Process:     return ostr<<"Signal Process    ";
  case btp::Hard_Decay:         return ostr<<"Hard Decay        ";
  case btp::Hard_Collision:     return ostr<<"Hard Collision    ";
  case btp::Soft_Collision:     return ostr<<"Soft Collision    "; 
  case btp::ME_PS_Interface_IS: return ostr<<"ME PS Interface   ";
  case btp::ME_PS_Interface_FS: return ostr<<"ME PS Interface   ";
  case btp::FS_Shower:          return ostr<<"FS Shower         ";
  case btp::IS_Shower:          return ostr<<"IS Shower         ";
  case btp::Beam:               return ostr<<"Beam              ";
  case btp::Bunch:              return ostr<<"Bunch             ";
  case btp::Fragmentation:      return ostr<<"Fragmentation     ";
  case btp::Cluster_Formation:  return ostr<<"Cluster Formation ";
  case btp::Cluster_Decay:      return ostr<<"Cluster Decay     ";
  case btp::Hadron_Decay:       return ostr<<"Hadron Decay      ";
  default:                      return ostr<<"Unknown           ";
  }
}

namespace ATOOLS {

  std::ostream& operator<<( std::ostream& ostr, const Blob & bl) {
#ifdef __GNUC__
#if __GNUC__ > 2
  std::ios_base::fmtflags flags=ostr.flags();
#else
  std::ios::fmtflags flags=ostr.flags();
#endif
#else
  std::ios::fmtflags flags=ostr.flags();
#endif
    ostr<<std::setw(4)<<std::setprecision(4);
    ostr<<"Blob ["<<bl.Status()<<"]( "<<bl.Id()<<", "<<bl.Type()<<", ";
    //    ostr<<"Blob ( "<<bl.Id()<<", "<<bl.Type()<<", ";
    if (bl.Beam() != -1) {
      ostr<<" from Beam "<<bl.Beam()<<", ";
    }
    ostr<<bl.NInP()<<" -> "<<bl.NOutP()<<" @ "<<bl.Position()<<std::endl;
    ostr<<"Incoming particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl.m_inparticles.begin();
	 part != bl.m_inparticles.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr<<"Outgoing particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl.m_outparticles.begin();
	 part != bl.m_outparticles.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr.setf(flags);
    return ostr;
  }


  std::ostream& operator<<( std::ostream& ostr,const  Blob * bl) {
#ifdef __GNUC__
#if __GNUC__ > 2
  std::ios_base::fmtflags flags=ostr.flags();
#else
  std::ios::fmtflags flags=ostr.flags();
#endif
#else
  std::ios::fmtflags flags=ostr.flags();
#endif
    ostr<<std::setw(4)<<std::setprecision(4);
    ostr<<"Blob ["<<bl->Status()<<"]( "<<bl->Id()<<", "<<bl->Type()<<", ";
    if (bl->Beam() != -1) {
      ostr<<" from Beam "<<bl->Beam()<<", ";
    }
    ostr<<bl->NInP()<<" -> "<<bl->NOutP()<<" @ "<<bl->Position()<<std::endl;
    ostr<<"Incoming particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl->m_inparticles.begin();
	 part != bl->m_inparticles.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr<<"Outgoing particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl->m_outparticles.begin();
	 part != bl->m_outparticles.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr.setf(flags);
    return ostr;
  }
}

Blob::Blob(const Vec4D _pos, const int _id) : 
  m_position(_pos), 
  m_id(_id), 
  m_weight(1.), 
  m_hasboost(false), 
  m_status(0), 
  m_beam(-1), 
  m_type(btp::Unspecified) {}

Blob::~Blob() {
  DeleteOwnedParticles();
  // delete data container
  ClearAllData();  
}

void Blob::AddToInParticles(Particle * Newp) {
  if (!Newp) return;
  m_inparticles.push_back( Newp );
  Newp->SetDecayBlob(this);
}

void Blob::AddToOutParticles(Particle * Newp) {
  if (!Newp) return;
  m_outparticles.push_back( Newp );
  Newp->SetProductionBlob(this);
}

Particle * Blob::InParticle(int _pos) {
  if (_pos>m_inparticles.size()-1 || _pos<0) { return NULL; }
  return m_inparticles[_pos];
}

Particle * Blob::OutParticle(int _pos) {
  if (_pos>m_outparticles.size()-1 || _pos<0) { return NULL; }
  return m_outparticles[_pos];
}

Particle * Blob::RemoveInParticle(int _pos,bool setit) {
  if (_pos>m_inparticles.size()-1 || _pos<0) { return NULL; }
  for (Particle_Vector::iterator part = m_inparticles.begin();
       part != m_inparticles.end(); ++part) {
    if ((*part)==m_inparticles[_pos]) {
      m_inparticles.erase(part);
      if (setit) (*part)->SetDecayBlob(NULL);
      return (*part);
    }
  }
  return NULL;
}

Particle * Blob::RemoveOutParticle(int _pos,bool setit) {
  if (_pos>m_outparticles.size()-1 || _pos<0) { return NULL; }
  for (Particle_Vector::iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    if ((*part)==m_outparticles[_pos]) {
      m_outparticles.erase(part);
      if (setit) (*part)->SetProductionBlob(NULL);
      return (*part);
    }
  }
  return NULL;
}


Particle * Blob::RemoveInParticle(Particle * _part,bool setit) {
  if (!_part) return 0;
  for (Particle_Vector::iterator part = m_inparticles.begin();
       part != m_inparticles.end(); ++part) {
    if ((*part)==_part) {
      Particle * p = (*part);
      m_inparticles.erase(part);
      if (setit) p->SetDecayBlob(NULL);
      return p;
    }
  }
  return NULL;
}

Particle * Blob::RemoveOutParticle(Particle * _part,bool setit) {
  if (!_part) return 0;
  for (Particle_Vector::iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    if ((*part)==_part) {
      Particle * p = (*part);
      m_outparticles.erase(part);
      if (setit) p->SetProductionBlob(NULL);
      return p;
    }
  }
  return NULL;
}


void Blob::DeleteInParticle(Particle * _part) {
  if (!_part) return;
  for (Particle_Vector::iterator part = m_inparticles.begin();
       part != m_inparticles.end(); ++part) {
    if ((*part)==_part) {
      m_inparticles.erase(part);
      if (_part->DecayBlob()==this) {
	if (_part->ProductionBlob()!=NULL) {
	  _part = _part->ProductionBlob()->RemoveOutParticle(_part);
	}
	delete _part;
	_part = NULL;
      }
      else msg.Out()<<"WARNING: particle not owned by the Blob asked to delete it"<<std::endl;
      return ;
    }
  }
}

void Blob::DeleteOutParticle(Particle * _part) {
  if (!_part) return;
  for (Particle_Vector::iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    if ((*part)==_part) {
      m_outparticles.erase(part);
      if (_part->ProductionBlob()==this) {
	if (_part->DecayBlob()!=NULL) {
	  _part = _part->DecayBlob()->RemoveInParticle(_part);
	}
	delete _part;
	_part = NULL;
      }
      else msg.Out()<<"WARNING: particle not owned by the Blob asked to delete it"<<std::endl;
      return ;
    }
  }
}



void Blob::DeleteOwnedParticles() {
  if (m_inparticles.empty() && m_outparticles.empty()) return;
  for (int i=m_inparticles.size()-1;i>=0;i--)  {
    DeleteInParticle(m_inparticles[i]);
  }
  for (int i=m_outparticles.size()-1;i>=0;i--) {
    DeleteOutParticle(m_outparticles[i]);
  }
  m_inparticles.clear();
  m_outparticles.clear();
}

Vec4D Blob::CheckMomentumConservation() {
  Vec4D sump = Vec4D(0.,0.,0.,0.);
  for (Particle_Vector::iterator part = m_inparticles.begin();
       part != m_inparticles.end(); ++part) {
    sump = sump + (*part)->Momentum();
  }
  for (Particle_Vector::iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    sump = sump + (-1.)*((*part)->Momentum());
  }
  return sump;
}

void Blob::BoostInCMS() {
  if (!m_hasboost) {
    Vec4D cm       = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<NInP();i++) cm = cm + InParticle(i)->Momentum();
    m_cms_boost = Poincare(cm);
    m_cms_vec   = cm;
  }
  for (int i=0;i<NInP();i++) 
    InParticle(i)->SetMomentum(m_cms_boost*InParticle(i)->Momentum());
  for (int i=0;i<NOutP();i++) 
    OutParticle(i)->SetMomentum(m_cms_boost*OutParticle(i)->Momentum());
  m_hasboost = true;
}

void Blob::BoostInLab() {
  if (!m_hasboost) {
    msg.Error()<<"Error in Blob::BoostInLab()."<<std::endl
	       <<"   Tried to boost back into unspecified system. Will just continue."<<std::endl;
  }
  Vec4D dummy;
  for (int i=0;i<NInP();i++) {
    dummy = InParticle(i)->Momentum();
    m_cms_boost.BoostBack(dummy);
    InParticle(i)->SetMomentum(dummy);
  }
  for (int i=0;i<NOutP();i++) { 
    dummy = OutParticle(i)->Momentum();
    m_cms_boost.BoostBack(dummy);
    OutParticle(i)->SetMomentum(dummy);
  }
}

void Blob::SetCMS() {
  m_cms_vec = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) m_cms_vec = m_cms_vec + InParticle(i)->Momentum();
}

void Blob::SetVecs() {
  m_cms_vec  = Vec4D(0.,0.,0.,0.);
  Vec4D  pos = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) {
    pos = pos + InParticle(i)->XDec();
  }
  for (int i=0;i<NOutP();i++) {
    m_cms_vec = m_cms_vec + OutParticle(i)->Momentum();
    pos = pos + OutParticle(i)->XProd();
  }
  m_position = 1./(NInP()+NOutP()) * pos;
}


void  Blob::AddData(const std::string name, Blob_Data_Base * data) 
{
  String_BlobDataBase_Map::iterator it=m_datacontainer.find(name);
  if (it==m_datacontainer.end()) {
    m_datacontainer[name]=data;
  }
  else {
    delete it->second;
    it->second=data;
  }
}

void Blob::ClearAllData() 
{
  for (String_BlobDataBase_Map::iterator it=m_datacontainer.begin();
       it!=m_datacontainer.end(); ++it) delete it->second;
  m_datacontainer.clear();
}


//=====================================================================



std::ostream& ATOOLS::operator<<( std::ostream& s, const Blob_Data_Base & bd) 
{
  bd>>s;
  return s;
}



Blob_Data_Base::~Blob_Data_Base()
{
}




template int Blob_Data_Base::Get<int>();
template long Blob_Data_Base::Get<long>();
template double Blob_Data_Base::Get<double>();
template std::string Blob_Data_Base::Get<std::string>();

template class Blob_Data<int>;
template class Blob_Data<long>;
template class Blob_Data<double>;
template class Blob_Data<std::string>;

