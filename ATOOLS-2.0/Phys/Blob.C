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

Parton * Blob::InParton(int n) {
  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    n--;
    if (n<0) return *part;
  }
  return 0;
}

Parton * Blob::OutParton(int n) {
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    n--;
    if (n<0) return *part;
  }
  return 0;
}

void Blob::AddToInPartons(Parton * Newp) {
  if (!Newp) return;
  m_inpartons.push_back( Newp );
  m_nin++;
}

void Blob::AddToOutPartons(Parton * Newp) {
  if (!Newp) return;
  m_outpartons.push_back( Newp );
  m_nout++;
}

Blob::~Blob() {
  DeleteOwnedPartons();
}

void Blob::DeleteOwnedPartons() {
  if ( (m_nin==0) && (m_nout==0) ) {
    if (m_inpartons.empty() && m_outpartons.empty()) {
      msg.Debugging()<<"Funny Partons in Blob !"<<std::endl; 
    }
    return;
  }
  if (m_inpartons.empty() && m_outpartons.empty()) {
    msg.Debugging()<<"Blob owns no more partons !"<<std::endl; 
    return;
  }
  for (Parton_Queue::iterator part = m_outpartons.begin();
       part != m_outpartons.end(); ++part) {
    delete *part;
    m_nout--;
  }
  m_outpartons.clear();

  for (Parton_Queue::iterator part = m_inpartons.begin();
       part != m_inpartons.end(); ++part) {
    delete *part;
    m_nin--;
  }
  m_inpartons.clear();
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


