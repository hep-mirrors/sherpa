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
      for (std::deque<Parton *>::const_iterator part = bl.InPartons.begin();
	   part != bl.InPartons.end(); ++part) {
	ostr<<*part<<std::endl;
      }
      //    }
    //    if (bl.NOutP() > 0) {
      ostr<<"Outgoing partons :"<<std::endl;
      for (std::deque<Parton *>::const_iterator part = bl.OutPartons.begin();
	   part != bl.OutPartons.end(); ++part) {
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
      for (std::deque<Parton *>::const_iterator part = bl->InPartons.begin();
	   part != bl->InPartons.end(); ++part) {
	ostr<<*part<<std::endl;
      }
      //    }
    //    if (bl->NOutP() > 0) {
      ostr<<"Outgoing partons :"<<std::endl;
      for (std::deque<Parton *>::const_iterator part = bl->OutPartons.begin();
	   part != bl->OutPartons.end(); ++part) {
	ostr<<*part<<std::endl;
      }
      //    }
    return ostr;
  }
}

Parton * Blob::InParton(int n) {
  for (std::deque<Parton *>::iterator part = InPartons.begin();
       part != InPartons.end(); ++part) {
    n--;
    if (n<0) return *part;
  }
  return 0;
}

Parton * Blob::OutParton(int n) {
  for (std::deque<Parton *>::iterator part = OutPartons.begin();
       part != OutPartons.end(); ++part) {
    n--;
    if (n<0) return *part;
  }
  return 0;
}

void Blob::AddToInPartons(Parton * Newp) {
  if (!Newp) return;
  InPartons.push_back( Newp );
  Nin++;
};

void Blob::AddToOutPartons(Parton * Newp) {
  if (!Newp) return;
  OutPartons.push_back( Newp );
  Nout++;
};

Blob::~Blob() {
  DeleteOwnedPartons();
};

void Blob::DeleteOwnedPartons() {
  if ( (Nin==0) && (Nout==0) ) {
    if (InPartons.empty() && OutPartons.empty()) {
      msg.Debugging()<<"Funny Partons in Blob !"<<std::endl; 
    }
    return;
  }
  if (InPartons.empty() && OutPartons.empty()) {
    msg.Debugging()<<"Blob owns no more partons !"<<std::endl; 
    return;
  }
  for (std::deque<Parton *>::iterator part = OutPartons.begin();
       part != OutPartons.end(); ++part) {
    delete *part;
    Nout--;
  }
  OutPartons.clear();

  for (std::deque<Parton *>::iterator part = InPartons.begin();
       part != InPartons.end(); ++part) {
    delete *part;
    Nin--;
  }
  InPartons.clear();
};

AMATOOLS::vec4d Blob::CheckMomentumConservation() {
  AMATOOLS::vec4d sump = AMATOOLS::vec4d(0.,0.,0.,0.);
  for (std::deque<Parton *>::iterator part = InPartons.begin();
       part != InPartons.end(); ++part) {
    sump = sump + (*part)->momentum();
  }
  for (std::deque<Parton *>::iterator part = OutPartons.begin();
       part != OutPartons.end(); ++part) {
    sump = sump + (-1.)*((*part)->momentum());
  }
  return sump;
};

void Blob::BoostInCMS() {
  AMATOOLS::vec4d cm       = vec4d(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) cm = cm + InParton(i)->momentum();
  cms_boost = AMATOOLS::Poincare(cm);

  cms                      = cm;
  for (int i=0;i<NInP();i++) 
    InParton(i)->set_momentum(cms_boost*InParton(i)->momentum());
  for (int i=0;i<NOutP();i++) 
    OutParton(i)->set_momentum( cms_boost*OutParton(i)->momentum());
};

void Blob::SetVecs() {
  cms  = vec4d(0.,0.,0.,0.);
  AMATOOLS::vec4d pos = vec4d(0.,0.,0.,0.);
  for (int i=0;i<NInP();i++) {
    pos = pos + InParton(i)->xdec();
  }
  for (int i=0;i<NOutP();i++) {
    cms = cms + OutParton(i)->momentum();
    pos = pos + OutParton(i)->xprod();
  }
  position = 1./(NInP()+NOutP()) * pos;
};

void Blob::BoostInLab() {
};

/*
Blob::Parton_Iterator::Parton_Iterator() {};

Blob::Parton_Iterator::Parton_Iterator(const Parton_Iterator& piter) {
  *this = piter;
};
*/

