#include "Blob.H"
#include "Parton.H"
#include "Poincare.H"
#include "Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ATOOLS {
  std::ostream& operator<<( std::ostream& ostr, const Blob & bl) {
    ostr<<std::setw(4)<<std::setprecision(4);
    ostr<<"Blob ( "<<bl.Id()<<", "<<bl.Type()<<", ";
    if (bl.Beam() != -1) {
      ostr<<" from Beam "<<bl.Beam()<<", ";
    }
    ostr<<bl.NInP()<<" -> "<<bl.NOutP()<<" @ "<<bl.Position()<<std::endl;
    ostr<<"Incoming partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl.m_inpartons.begin();
	 part != bl.m_inpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr<<"Outgoing partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl.m_outpartons.begin();
	 part != bl.m_outpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    return ostr;
  }


  std::ostream& operator<<( std::ostream& ostr,const  Blob * bl) {
    ostr<<std::setw(4)<<std::setprecision(4);
    ostr<<"Blob ( "<<bl->Id()<<", "<<bl->Type()<<", ";
    if (bl->Beam() != -1) {
      ostr<<" from Beam "<<bl->Beam()<<", ";
    }
    ostr<<bl->NInP()<<" -> "<<bl->NOutP()<<" @ "<<bl->Position()<<std::endl;
    ostr<<"Incoming partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl->m_inpartons.begin();
	 part != bl->m_inpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    ostr<<"Outgoing partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl->m_outpartons.begin();
	 part != bl->m_outpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    return ostr;
  }
}

Blob::~Blob() {
  DeleteOwnedPartons();
}

void Blob::AddToInPartons(Parton * Newp) {
  if (!Newp) return;
  m_inpartons.push_back( Newp );
}

void Blob::AddToOutPartons(Parton * Newp) {
  if (!Newp) return;
  m_outpartons.push_back( Newp );
}

Parton * Blob::InParton(int _pos) {
  if (_pos>m_inpartons.size()-1 || _pos<0) { return NULL; }
  return m_inpartons[_pos];
}

Parton * Blob::OutParton(int _pos) {
  if (_pos>m_outpartons.size()-1 || _pos<0) { return NULL; }
  return m_outpartons[_pos];
}

Parton * Blob::RemoveInParton(int _pos) {
  if (_pos>m_inpartons.size()-1 || _pos<0) { return NULL; }
  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    if ((*part)==m_inpartons[_pos]) {
      m_inpartons.erase(part);
      (*part)->SetDecayBlob(NULL);
      return (*part);
    }
  }
  return NULL;
}

Parton * Blob::RemoveOutParton(int _pos) {
  if (_pos>m_outpartons.size()-1 || _pos<0) { return NULL; }
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    if ((*part)==m_outpartons[_pos]) {
      m_outpartons.erase(part);
      (*part)->SetProductionBlob(NULL);
      return (*part);
    }
  }
  return NULL;
}


Parton * Blob::RemoveInParton(Parton * _part) {
  if (!_part) return 0;
  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    if ((*part)==_part) {
      Parton * p = (*part);
      m_inpartons.erase(part);
      p->SetDecayBlob(NULL);
      return p;
    }
  }
  return NULL;
}

Parton * Blob::RemoveOutParton(Parton * _part) {
  if (!_part) return 0;
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    if ((*part)==_part) {
      Parton * p = (*part);
      m_outpartons.erase(part);
      p->SetProductionBlob(NULL);
      return p;
    }
  }
  return NULL;
}


void Blob::DeleteInParton(Parton * _part) {
  if (!_part) return;
  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    if ((*part)==_part) {
      m_inpartons.erase(part);
      if (_part->DecayBlob()==this) {
	if (_part->ProductionBlob()!=NULL) {
	  _part = _part->ProductionBlob()->RemoveOutParton(_part);
	}
	delete _part;
	_part = NULL;
      }
      // else msg.Out()<<"WARNING: parton not owned by the Blob asked to delete it"<<std::endl;
      return ;
    }
  }
}

void Blob::DeleteOutParton(Parton * _part) {
  if (!_part) return;
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    if ((*part)==_part) {
      m_outpartons.erase(part);
      if (_part->ProductionBlob()==this) {
	if (_part->DecayBlob()!=NULL) {
	  _part = _part->DecayBlob()->RemoveInParton(_part);
	}
	delete _part;
	_part = NULL;
      }
      // else msg.Out()<<"WARNING: parton not owned by the Blob ask to delete it"<<std::endl;
      return ;
    }
  }
}



void Blob::DeleteOwnedPartons() {
  if (m_inpartons.empty() && m_outpartons.empty()) return;
  for (int i=m_inpartons.size()-1;i>=0;i--)  {
    DeleteInParton(m_inpartons[i]);
  }
  for (int i=m_outpartons.size()-1;i>=0;i--) {
    DeleteOutParton(m_outpartons[i]);
  }
  m_inpartons.clear();
  m_outpartons.clear();
}

Vec4D Blob::CheckMomentumConservation() {
  Vec4D sump = Vec4D(0.,0.,0.,0.);
  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    sump = sump + (*part)->Momentum();
  }
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    sump = sump + (-1.)*((*part)->Momentum());
  }
  return sump;
}

void Blob::BoostInCMS() {
  Vec4D cm       = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) cm = cm + InParton(i)->Momentum();
  m_cms_boost = Poincare(cm);

  m_cms_vec                      = cm;
  for (int i=0;i<NInP();i++) 
    InParton(i)->SetMomentum(m_cms_boost*InParton(i)->Momentum());
  for (int i=0;i<NOutP();i++) 
    OutParton(i)->SetMomentum(m_cms_boost*OutParton(i)->Momentum());
}

void Blob::SetCMS() {
  m_cms_vec = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) m_cms_vec = m_cms_vec + InParton(i)->Momentum();
}

void Blob::SetVecs() {
  m_cms_vec  = Vec4D(0.,0.,0.,0.);
  Vec4D  pos = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) {
    pos = pos + InParton(i)->XDec();
  }
  for (int i=0;i<NOutP();i++) {
    m_cms_vec = m_cms_vec + OutParton(i)->Momentum();
    pos = pos + OutParton(i)->XProd();
  }
  m_position = 1./(NInP()+NOutP()) * pos;
}

void Blob::BoostInLab() {

}


