#include "Blob.H"
#include "Parton.H"
#include "Poincare.H"
#include "Message.H"


using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;

namespace APHYTOOLS {
  std::ostream& operator<<( std::ostream& ostr, const Blob & bl) {
    ostr<<"Blob ( "<<bl.Id()<<", "<<bl.Type()<<" ), ";
    if (bl.Beam() != -1) {
      ostr<<" from Beam "<<bl.Beam()<<", ";
    }
    ostr<<std::endl<<" #in = "<<bl.NInP()<<", #out = "<<bl.NOutP()<<std::endl
	<<" at "<<bl.Position()<<std::endl<<" with "<<bl.CMS()<<std::endl;
    //    if (bl.NInP() > 0) {
      ostr<<"Incoming partons :"<<std::endl;
      for (Parton_Queue::const_iterator part = bl.m_inpartons.begin();
	   part != bl.m_inpartons.end(); ++part) {
	ostr<<*part<<std::endl;
      }
      //    }
    //    if (bl.NOutP() > 0) {
      ostr<<"Outgoing partons :"<<std::endl;
      for (Parton_Queue::const_iterator part = bl.m_outpartons.begin();
	   part != bl.m_outpartons.end(); ++part) {
	ostr<<*part<<std::endl;
      }
      //    }
    return ostr;
  }


  std::ostream& operator<<( std::ostream& ostr,const  Blob * bl) {
    ostr<<"Blob ( "<<bl->Id()<<", "<<bl->Type()<<" ), ";
    if (bl->Beam() != -1) {
      ostr<<" from Beam "<<bl->Beam()<<", ";
    }
    ostr<<std::endl<<" #in = "<<bl->NInP()<<", #out = "<<bl->NOutP()<<std::endl
	<<" at "<<bl->Position()<<std::endl
	<<" with "<<bl->CMS()<<std::endl;
    //    if (bl->NInP() > 0) {
    ostr<<"Incoming partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl->m_inpartons.begin();
	 part != bl->m_inpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    //    }
    //    if (bl->NOutP() > 0) {
    ostr<<"Outgoing partons :"<<std::endl;
    for (Parton_Queue::const_iterator part = bl->m_outpartons.begin();
	 part != bl->m_outpartons.end(); ++part) {
      ostr<<*part<<std::endl;
    }
    //    }
    return ostr;
  }
}

Blob::~Blob() {
  DeleteOwnedPartons();
  std::cout<<"Deleted all partons in blob : "<<m_id<<std::endl;
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
      m_inpartons.erase(part);
      (*part)->SetDecayBlob(NULL);
      return (*part);
    }
  }
  return NULL;
}

Parton * Blob::RemoveOutParton(Parton * _part) {
  if (!_part) return 0;
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    if ((*part)==_part) {
      m_outpartons.erase(part);
      (*part)->SetProductionBlob(NULL);
      return (*part);
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
      if (_part->ProductionBlob()!=NULL) {
	_part = _part->ProductionBlob()->RemoveOutParton(_part);
      }
      delete _part;
      _part = NULL;
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
      if (_part->DecayBlob()!=NULL) {
	_part = _part->DecayBlob()->RemoveInParton(_part);
      }
      delete _part;
      _part = NULL;
      return ;
    }
  }
}



void Blob::DeleteOwnedPartons() {
  if (m_inpartons.empty() && m_outpartons.empty()) return;
  std::cout<<"In DeleteOwnedPartons() for blob no. : "<<m_id<<" with "
      <<m_inpartons.size()<<" -> "<<m_outpartons.size()<<std::endl;
  for (int i=m_inpartons.size()-1;i>=0;i--)  {
    std::cout<<"Try to delete inparton : "<<i<<"/"<<m_inpartons.size()-1<<std::endl;
    DeleteInParton(m_inpartons[i]);
    std::cout<<"Succeeded."<<std::endl;
  }
  for (int i=m_outpartons.size()-1;i>=0;i--) {
    std::cout<<"Try to delete outparton : "<<i<<"/"<<m_outpartons.size()-1<<std::endl;
    DeleteOutParton(m_outpartons[i]);
    std::cout<<"Succeeded."<<std::endl;
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


