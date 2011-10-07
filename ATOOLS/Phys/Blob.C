#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

std::ostream& ATOOLS::operator<<(std::ostream& ostr, const btp::code btpc) {
  switch (btpc) {
  case btp::Unspecified:        return ostr<<"Unspecified       ";
  case btp::Signal_Process:     return ostr<<"Signal Process    ";
  case btp::Hard_Decay:         return ostr<<"Hard Decay        ";
  case btp::Hard_Collision:     return ostr<<"Hard Collision    ";
  case btp::Soft_Collision:     return ostr<<"Soft Collision    "; 
  case btp::Shower:             return ostr<<"Shower            ";
  case btp::QED_Radiation:      return ostr<<"QED Radiation     ";
  case btp::Beam:               return ostr<<"Beam              ";
  case btp::Bunch:              return ostr<<"Bunch             ";
  case btp::Fragmentation:      return ostr<<"Fragmentation     ";
  case btp::Cluster_Formation:  return ostr<<"Cluster Formation ";
  case btp::Cluster_Decay:      return ostr<<"Cluster Decay     ";
  case btp::Hadron_Decay:       return ostr<<"Hadron Decay      ";
  case btp::Hadron_Mixing:      return ostr<<"Hadron Mixing     ";
  case btp::Hadron_To_Parton:   return ostr<<"Hadron-To-Partons ";
  default:                      return ostr<<"Unknown           ";
  }
}

namespace ATOOLS {
  int Blob::s_totalnumber=0;
  long unsigned int Blob::s_currentnumber=0;
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
    ostr<<bl.NInP()<<" -> "<<bl.NOutP()<<" @ "<<bl.Position()<<std::endl;
    ostr<<"Incoming particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl.m_inparticles.begin();
	 part != bl.m_inparticles.end(); ++part) {
      ostr<<**part<<std::endl;
    }
    ostr<<"Outgoing particles :"<<std::endl;
    for (Particle_Vector::const_iterator part = bl.m_outparticles.begin();
	 part != bl.m_outparticles.end(); ++part) {
      ostr<<**part<<std::endl;
    }
    if (bl.m_datacontainer.size()>0) {
      ostr<<"Data_Container:"<<std::endl;
      for (String_BlobDataBase_Map::const_iterator it=bl.m_datacontainer.begin();
	   it!=bl.m_datacontainer.end(); ++it) {
	ostr<<"   * "<<it->first<<" ("<<*(it->second)<<")"<<std::endl;
      }
    }
    ostr.setf(flags);
    return ostr;
  }

}

Blob::Blob(const Vec4D _pos, const int _id) : 
  m_position(_pos), m_id(_id), m_weight(1.), m_status(blob_status::inactive), 
  m_hasboost(false), 
  m_type(btp::Unspecified), m_typespec(std::string("none")) 
{ ++s_totalnumber; }

Blob::Blob(const Blob * blob,const bool copyparts) :
  m_position(blob->m_position), m_id(blob->m_id), m_weight(blob->m_weight),
  m_status(blob->m_status),
  m_type(blob->m_type), m_typespec(blob->m_typespec),
  m_cms_vec(blob->m_cms_vec), m_cms_boost(Poincare(m_cms_vec))
{
  ++s_totalnumber;
  if (copyparts) {
  for (int i=0;i<blob->NInP();i++)  AddToInParticles(new Particle((*blob->ConstInParticle(i))));
  Particle * part(NULL);
  for (int i=0;i<blob->NOutP();i++) {
    part = new Particle((*blob->ConstOutParticle(i)));
    part->SetStatus(part_status::active);
    AddToOutParticles(part);
  }
  }
  for (String_BlobDataBase_Map::const_iterator it=blob->m_datacontainer.begin(); it!=blob->m_datacontainer.end(); ++it) {
    AddData(it->first, it->second->ClonePtr());
  }
}

Blob::~Blob() {
  DeleteOwnedParticles();
  // delete data container
  ClearAllData();  
  --s_totalnumber;
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
  if (_pos>(int)m_inparticles.size()-1 || _pos<0) { return NULL; }
  return m_inparticles[_pos];
}

Particle * Blob::OutParticle(int _pos) {
  if (_pos>(int)m_outparticles.size()-1 || _pos<0) { return NULL; }
  return m_outparticles[_pos];
}

const Particle *Blob::ConstInParticle(const size_t i) const
{
  if (i>m_inparticles.size()-1 || i<0) return NULL;
  return m_inparticles[i];
}

const Particle *Blob::ConstOutParticle(const size_t i) const 
{
  if (i>m_outparticles.size()-1 || i<0) return NULL; 
  return m_outparticles[i];
}

Particle * Blob::RemoveInParticle(int _pos,bool setit) {
  if (_pos>(int)m_inparticles.size()-1 || _pos<0) { return NULL; }
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
  if (_pos>(int)m_outparticles.size()-1 || _pos<0) { return NULL; }
  for (Particle_Vector::iterator pit = m_outparticles.begin();
       pit != m_outparticles.end(); ++pit) {
    Particle * part(*pit);
    if (part==m_outparticles[_pos]) {
      m_outparticles.erase(pit);
      if (setit) part->SetProductionBlob(NULL);
      return part;
    }
  }
  return NULL;
}

void Blob::RemoveInParticles(const int all) 
{
  for (Particle_Vector::iterator part=m_inparticles.begin();
       part!=m_inparticles.end();) {
    if ((all==-1&&(*part)->ProductionBlob()==NULL) ||
	all==0 ||
	(all==1&&(*part)->ProductionBlob()!=NULL)) {
      (*part)->SetDecayBlob(NULL);
      part=m_inparticles.erase(part);
    }
    else {
      ++part;
    }
  }
}

void Blob::RemoveOutParticles(const int all) 
{
  for (Particle_Vector::iterator part=m_outparticles.begin();
       part!=m_outparticles.end();) {
    if ((all==-1&&(*part)->DecayBlob()==NULL) ||
	all==0 ||
	(all==1&&(*part)->DecayBlob()!=NULL)) {
      (*part)->SetProductionBlob(NULL);
      part=m_outparticles.erase(part);
    }
    else {
      ++part;
    }
  }
}

void Blob::DeleteInParticles(const int all) 
{
  for (Particle_Vector::iterator part=m_inparticles.begin();
       part!=m_inparticles.end();) {
    if ((all==-1&&(*part)->ProductionBlob()==NULL) ||
	all==0 ||
	(all==1&&(*part)->ProductionBlob()!=NULL)) {
      if ((*part)->ProductionBlob()!=NULL) 
	(*part)->ProductionBlob()->RemoveOutParticle(*part);
      (*part)->SetDecayBlob(NULL);
      delete *part;
      part=m_inparticles.erase(part);
    }
    else {
      ++part;
    }
  }
}

void Blob::DeleteOutParticles(const int all) 
{
  for (Particle_Vector::iterator part=m_outparticles.begin();
       part!=m_outparticles.end();) {
    if ((all==-1&&(*part)->DecayBlob()==NULL) ||
	all==0 ||
	(all==1&&(*part)->DecayBlob()!=NULL)) {
      if ((*part)->DecayBlob()!=NULL) 
	(*part)->DecayBlob()->RemoveInParticle(*part);
      (*part)->SetProductionBlob(NULL);
      delete *part;
      part=m_outparticles.erase(part);
    }
    else {
      ++part;
    }
  }
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
      if (_part->DecayBlob()==this) {
	if (_part->ProductionBlob()!=NULL) _part->ProductionBlob()->RemoveOutParticle(_part);
	delete _part;
	_part = NULL;
      }
      else {
	msg_Out()<<"WARNING in "<<METHOD<<":"<<std::endl
		 <<"   particle not owned by the Blob asked to delete it"<<std::endl
		 <<"   "<<(*_part)<<std::endl;
      }
      m_inparticles.erase(part);
      return ;
    }
  }
}

void Blob::RemoveOwnedParticles(const bool del)
{
  for (size_t i=0;i<m_inparticles.size();++i) {
    if (m_inparticles[i]->ProductionBlob()==NULL) {
      if (del) delete m_inparticles[i];
    }
    else m_inparticles[i]->SetDecayBlob(NULL);
  }
  m_inparticles.clear();
  for (size_t i=0;i<m_outparticles.size();++i) {
    if (m_outparticles[i]->DecayBlob()==NULL) {
      if (del) delete m_outparticles[i];
    }
    else m_outparticles[i]->SetProductionBlob(NULL);
  }
  m_outparticles.clear();
}

void Blob::DeleteOutParticle(Particle * _part) {
  if (!_part) return;
  for (Particle_Vector::iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    if ((*part)==_part) {
      m_outparticles.erase(part);
      if (_part->ProductionBlob()==this) {
	if (_part->DecayBlob()!=NULL) _part->DecayBlob()->RemoveInParticle(_part);
	delete _part;
	_part = NULL;
      }
      else {
	msg_Out()<<"WARNING in "<<METHOD<<":"<<std::endl
		 <<"   particle not owned by the Blob asked to delete it"<<std::endl
		 <<"   "<<(*_part)<<std::endl;
      }
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

Vec4D Blob::CheckMomentumConservation() const {
  Vec4D sump = Vec4D(0.,0.,0.,0.);
  for (Particle_Vector::const_iterator part = m_inparticles.begin();
       part != m_inparticles.end(); ++part) {
    if ((*part)->Info()=='F' && m_type==btp::Shower) sump = sump + (-1.) * (*part)->Momentum();
                                                else sump = sump + (*part)->Momentum();
  }
  for (Particle_Vector::const_iterator part = m_outparticles.begin();
       part != m_outparticles.end(); ++part) {
    if ((*part)->Info()=='I' && m_type==btp::Shower) sump = sump + (*part)->Momentum();
                         else sump = sump + (-1.)*((*part)->Momentum());
  }
  return sump;
}

bool Blob::MomentumConserved() {
  Vec4D cms_vec = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) cms_vec = cms_vec + InParticle(i)->Momentum();
  Vec4D mc = CheckMomentumConservation();
  double accu=1e-6*cms_vec[0];
  if(dabs(mc[0])>accu||dabs(mc[1])>accu||dabs(mc[2])>accu||dabs(mc[3])>accu) {
    return false;
  }
  return true;
}

void Blob::Boost(const Poincare& boost) {
  for (int i=0;i<NInP();i++)
    InParticle(i)->SetMomentum(boost*InParticle(i)->Momentum());
  for (int i=0;i<NOutP();i++)
    OutParticle(i)->SetMomentum(boost*OutParticle(i)->Momentum());
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
    msg_Error()<<"Error in Blob::BoostInLab()."<<std::endl
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


void  Blob::SetId(const int _id) { 
  if (_id<0) m_id = -_id;
        else m_id = ++s_currentnumber; 
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

long int Blob_Data_Base::s_number(0);

Blob_Data_Base::Blob_Data_Base()
{
  ++s_number;
}

Blob_Data_Base::Blob_Data_Base(const Blob_Data_Base &base)
{
  ++s_number;
}

Blob_Data_Base::~Blob_Data_Base()
{
  --s_number;
}

template <class Type>
Blob_Data<Type>::~Blob_Data() 
{
}


template int &Blob_Data_Base::Get<int>();
template size_t &Blob_Data_Base::Get<size_t>();
template long &Blob_Data_Base::Get<long>();
template double &Blob_Data_Base::Get<double>();
template std::string &Blob_Data_Base::Get<std::string>();
template std::vector<double> &Blob_Data_Base::Get<std::vector<double> >();
template std::vector<int> &Blob_Data_Base::Get<std::vector<int> >();
template Vec4D &Blob_Data_Base::Get<Vec4D>();

template void Blob_Data_Base::Set<int>(const int &data);
template void Blob_Data_Base::Set<size_t>(const size_t &data);
template void Blob_Data_Base::Set<long>(const long &data);
template void Blob_Data_Base::Set<double>(const double &data);
template void Blob_Data_Base::Set<std::string>(const std::string &data);
template void Blob_Data_Base::Set<std::vector<double> >(const std::vector<double> &data);
template void Blob_Data_Base::Set<std::vector<int> >(const std::vector<int> &data);
template void Blob_Data_Base::Set<Vec4D>(const Vec4D &data);

template class Blob_Data<int>;
template class Blob_Data<size_t>;
template class Blob_Data<long>;
template class Blob_Data<double>;
template class Blob_Data<std::string>;
template class Blob_Data<std::vector<double> >;
template class Blob_Data<std::vector<int> >;
template class Blob_Data<Vec4D>;

void Blob::SwapInParticles(const size_t i, const size_t j) 
{
  if (i<m_inparticles.size() && j<m_inparticles.size()) {
    ATOOLS::Particle *help=m_inparticles[j];
    m_inparticles[j]=m_inparticles[i];
    m_inparticles[i]=help;
  }
}

void Blob::SwapOutParticles(const size_t i, const size_t j) 
{
  if (i<m_outparticles.size() && j<m_outparticles.size()) {
    ATOOLS::Particle *help=m_outparticles[j];
    m_outparticles[j]=m_outparticles[i];
    m_outparticles[i]=help;
  }
}


bool Blob::IsConnectedTo(const btp::code &type,
			 std::set<const Blob*> &checked) const
{
  if (this==NULL || checked.find(this)!=checked.end()) return false;
  checked.insert(this);
  if (Type()==type) return true;
  for (int i(0);i<NOutP();++i) 
    if (ConstOutParticle(i)->DecayBlob()->IsConnectedTo(type,checked)) 
      return true;
  for (int i(0);i<NInP();++i) 
    if (ConstInParticle(i)->ProductionBlob()->IsConnectedTo(type,checked)) 
      return true;
  return false;
}

bool Blob::IsConnectedTo(const btp::code &type) const
{
  std::set<const Blob*> checked;
  return IsConnectedTo(type,checked);
}

Blob* Blob::UpstreamBlob() const
{
  if(NInP()==0) return NULL;
  Blob* upstream = ConstInParticle(0)->ProductionBlob();
  for (int i(1);i<NInP();++i) {
    if(ConstInParticle(i)->ProductionBlob()!=upstream) return NULL;
  }
  return upstream;
}

Blob* Blob::DownstreamBlob() const
{
  if(NOutP()==0) return NULL;
  Blob* downstream = ConstOutParticle(0)->DecayBlob();
  for (int i(1);i<NOutP();++i) {
    if(ConstOutParticle(i)->DecayBlob()!=downstream) return NULL;
  }
  return downstream;
}
