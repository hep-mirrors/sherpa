//bof
//Version: 2 ADICIC++-0.0/2005/01/31

//Implementation of Chain.H.



#ifdef __GNUC__
#if __GNUC__ >2
#include <ios>
#endif
#endif
#include <iomanip>
#include <string>
#include <sstream>
#include "Message.H"
#include "Chain.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



ostream& ADICIC::operator<<(ostream& ost, const ADICIC::Chain& cha) {
  unsigned i=0;
  ost<<"\n"<<om::bold<<string(114,'=')<<" "<<cha.Name<<" | "
     <<cha.Memo<<" =="<<om::reset<<"\n";
  for(list<Dipole*>::const_iterator it=cha.DipolePointerList().begin();
      it!=cha.DipolePointerList().end(); ++it) {
    ++i;
    if(*it==&cha.ChainRoot()) ost<<om::bold<<om::red;
    else ost<<om::bold<<om::blue;
    ost<<setiosflags(ios::left)<<setw(4)<<i<<om::reset
       <<resetiosflags(ios::left)<<(**it)<<"\n";
  }
  ost<<om::bold<<" "<<string(118,'-')<<om::reset<<"\n";
  ost<<"Chain type: ";
  switch(cha.ChainType()) {
  case 1 : ost<<"ring"; break;
  case 0 : ost<<"line"; break;
  default: ost<<"incorrect";
  }
  ost<<"    Chain state: ";
  switch(cha.Status()) {
  case -1: ost<<"blocked"; break;
  case 0 : ost<<"off"; break;
  default: ost<<"on";
  }
  ost<<"    Chain handler: "<<cha.IsHandled()<<"\n";
  ost<<setiosflags(ios::left)
     <<"Scale="<<setw(12)<<cha.LastScale()
     <<"Mass="<<setw(12)<<cha.Mass()
     <<"SqrMass="<<setw(12)<<cha.InvMass()
     <<"P="<<cha.Momentum();
  ost<<resetiosflags(ios::left)<<"\n";
  if(Chain::s_print) {
    ost<<om::bold<<" "<<string(118,'-')<<om::reset<<"\n";
    if(cha.BranchPointer()) ost<<om::bold<<"Branch:"<<om::reset<<"\n"
			       <<cha.BranchPointer()->Parton<<"\n";
    if(cha.AntibranchPointer()) ost<<om::bold<<"Antibranch:"<<om::reset<<"\n"
				   <<cha.AntibranchPointer()->Parton<<"\n";
    if(!cha.GlubranchPointerList().empty())
      ost<<om::bold<<"Glubranches:"<<om::reset<<"\n";
    for(list<Dipole::Glubranch*>::const_iterator
	  gut=cha.GlubranchPointerList().begin();
	gut!=cha.GlubranchPointerList().end(); ++gut) {
      if(*gut==cha.FirstGlubranchPointer())
	ost<<(*gut)->Parton<<"  "<<om::red<<"g_root"<<om::reset<<"\n";
      else
	ost<<(*gut)->Parton<<"\n";
    }
    ost<<om::bold<<" "<<string(118,'-')<<om::reset<<"\n";
    ost<<"Number of branches = "<<cha.ParticleNumber()<<"\n";
    ost<<"Number of  dipoles = "<<cha.DipoleNumber()<<"\n";
  }
  ost<<om::bold<<string(114,'=')<<" "<<cha.Name<<" | "
     <<cha.Memo<<" =="<<om::reset;
  return ost;
}



//=============================================================================



const Chain::Initiator::Simple_EpEm Chain::Initiator::simple_epem={};



//=============================================================================



bool Chain::s_print=false;
bool& Chain::MoreCout=Chain::s_print;

size_t Chain::s_count=0;
size_t Chain::s_maxcount=0;
const size_t& Chain::InStore=Chain::s_count;
const size_t& Chain::MaxCount=Chain::s_maxcount;



//=============================================================================



Chain::Chain()
  : m_name(++s_maxcount), m_memo(0), Name(m_name), Memo(m_memo) {
  ++s_count;
  varset.Init();    //This is a clear chain.
}





Chain::Chain(const Chain& cha)
  : m_name(++s_maxcount), m_memo(0), Name(m_name), Memo(m_memo) {
  ++s_count;
  varset.Init();
  //Later: Additional Chain initialization methods will cause further
  //modifications here!
  Type type=cha.ChainType();
  if(type==incorrect) return;//////////////////////////////////////////////////
  //Chain is not empty here!
  this->Copy(cha,type);
  m_memo=m_name;
}





Chain::Chain(const Dipole::Branch& ban, const Dipole::Antibranch& ati,
	     const Initiator::Simple_EpEm)
  : m_name(++s_maxcount), m_memo(m_name), Name(m_name), Memo(m_memo) {

  ++s_count;

  varset.f_clear=false;
  varset.p_hdl=NULL;
  varset.p_1glu=NULL;
  varset.l_glub.clear();
  varset.l_dip.clear();

  varset.p_quab=new Dipole::Branch(ban);    //Copying!!
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);
  varset.p_root=root;

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

}





Chain::Chain(const Dipole::Glubranch& glut, const Dipole::Glubranch& glub,
	     const Initiator::Simple_EpEm)
  : m_name(++s_maxcount), m_memo(m_name), Name(m_name), Memo(m_memo) {

  ++s_count;

  if(&glut==&glub) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Gluon Chain initiation from a single gluon!\n"<<endl;
    //exit(1);
  }

  varset.f_clear=false;
  varset.p_hdl=NULL;
  varset.p_quab=NULL;
  varset.p_atib=NULL;
  varset.l_glub.clear();
  varset.l_dip.clear();

  Dipole::Glubranch* top=new Dipole::Glubranch(glut);
  assert(top);
  varset.l_glub.push_front(top);
  varset.p_1glu=top;
  Dipole::Glubranch* bot=new Dipole::Glubranch(glub);
  assert(bot);
  varset.l_glub.push_back(bot);

  //First Dipole.
  Dipole* root=new Dipole(*top,*bot);
  assert(root);
  varset.l_dip.push_front(root);
  varset.p_root=root;
  //Second Dipole.
  root=new Dipole(*bot,*top);
  assert(root);
  varset.l_dip.push_back(root);

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

}





Chain::~Chain() {

  if(varset.p_hdl) {    //Check for Chain_Handler.
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Detaching Chain_Handler from Chain!\n"<<endl;
    if(varset.p_hdl->IsDockedAt(*this)==false) {
      cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
      assert(varset.p_hdl->IsDockedAt(*this));
    }
    Chain_Handler* phand=varset.p_hdl;
    varset.p_hdl=NULL;
    phand->DetachChain(this);
  }

  for(list<Dipole*>::iterator dit=varset.l_dip.begin();
      dit!=varset.l_dip.end(); ++dit)
    if(*dit) delete (*dit);
  for(list<Dipole::Glubranch*>::iterator git=varset.l_glub.begin();
      git!=varset.l_glub.end(); ++git)
    if(*git) delete (*git);
  if(varset.p_quab) delete varset.p_quab;
  if(varset.p_atib) delete varset.p_atib;

  varset.l_dip.clear();
  varset.l_glub.clear();

  --s_count;
  if(m_name && m_name==s_maxcount) --s_maxcount;
  if(!s_count) s_maxcount=0;

}



//=============================================================================



Chain& Chain::operator=(const Chain& cha) {
  if(this==&cha) return *this;
  //Later: Additional Chain initialization methods will cause further
  //modifications here!
  Type type=cha.ChainType();
  if(type==incorrect) return *this;////////////////////////////////////////////
  if(this->Clear()==false) return *this;
  this->Copy(cha,type);
  m_memo=m_name;
  return *this;
}





//const bool Chain::operator==(const Chain& cha) const {
//  cerr<<"\nSorry :o( Method has not been implemented yet.\n";
//  assert(0); exit(1);
//}



//=============================================================================



const Chain::Type Chain::ChainType() const {
  if(varset.l_dip.empty()) return Chain::incorrect;
  if(varset.l_dip.front()->GetTopBranchPointer().operator->()==varset.p_quab
     &&
     varset.l_dip.back()->GetBotBranchPointer().operator->()==varset.p_atib)
    return Chain::line;
  if(varset.l_dip.front()->GetTopBranchPointer().operator->()==varset.p_1glu
     &&
     varset.l_dip.back()->GetBotBranchPointer().operator->()==varset.p_1glu)
    return Chain::ring;
  return Chain::incorrect;
}





const bool Chain::CheckMomentumConservation(ATOOLS::Vec4D& sum) const {
  sum=Vec4D();
  if(varset.p_quab) sum+=varset.p_quab->Momentum();
  if(varset.p_atib) sum+=varset.p_atib->Momentum();
  for(list<Dipole::Glubranch*>::const_iterator gut=varset.l_glub.begin();
      gut!=varset.l_glub.end(); ++gut)
    sum+=(*gut)->Momentum();
  if(sum==varset.m_momentum) return true;
  return false;
}





const bool Chain::ExtractPartons(Particle_List& parlist) const {

  //The chain entries are appended to the parlist.

  Type type=this->ChainType();
  if(type==incorrect) return false;

  list<Dipole*>::const_iterator dit=varset.l_dip.begin();
  Particle* par=new Particle( (*dit)->GetTopBranchPointer()->Parton );
  parlist.push_back(par);
  Particle_List::iterator pit=parlist.end(); --pit;
  ++dit;
  for(; dit!=varset.l_dip.end(); ++dit) {
    par=new Particle( (*dit)->GetTopBranchPointer()->Parton );
    parlist.push_back(par);
  }

  if(type==line) {
    --dit;
    Particle* par=new Particle( (*dit)->GetBotBranchPointer()->Parton );
    parlist.push_back(par);

    const size_t num=parlist.size(); assert(num>1);
    if(parlist.size()>2) {
      unsigned upflow=(*pit)->GetFlow(1);
      if(upflow==0) { (*pit)->SetFlow(1,-1); upflow=(*pit)->GetFlow(1);}
      ++pit;
      Particle_List::iterator pot=parlist.end();
      --pot;
      for(; pit!=pot; ++pit) {
	(*pit)->SetFlow(2,upflow);
	(*pit)->SetFlow(1,-1);
	upflow=(*pit)->GetFlow(1);
      }
      (*pit)->SetFlow(2,upflow);
    } else {
      unsigned upflow=(*pit)->GetFlow(1);
      if(upflow==0) { (*pit)->SetFlow(1,-1); upflow=(*pit)->GetFlow(1);}
      ++pit;
      if((*pit)->GetFlow(2)!=upflow) (*pit)->SetFlow(2,upflow);
    }

#ifdef CHAIN_OUTPUT
    cout<<parlist<<endl;
#endif

    return true;
  }

  //else ... type==ring ... No flow setting so far.

#ifdef CHAIN_OUTPUT    //For testing, it is suitable to comment it out.
  cout<<parlist<<endl;
#endif    //For testing, it is suitable to comment it out.

  return true;

}





const bool Chain::Clear() {
  if(varset.p_hdl) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Clearing a Chain "
	<<"manipulated by a Chain_Handler is not permitted!\n"<<endl;
    return false;
  }
  for(list<Dipole*>::iterator dit=varset.l_dip.begin();
      dit!=varset.l_dip.end(); ++dit) {
    if(*dit) delete (*dit);
  }
  for(list<Dipole::Glubranch*>::iterator git=varset.l_glub.begin();
      git!=varset.l_glub.end(); ++git) {
    if(*git) delete (*git);
  }
  if(varset.p_quab) delete varset.p_quab;
  if(varset.p_atib) delete varset.p_atib;
  varset.Init();    //This is a clear chain right now.
  m_memo=0;
  return true;
}





const bool Chain::Initialize(const Dipole::Branch& ban,
			     const Dipole::Antibranch& ati) {

  //if(!IsClear()) return false;
  if(varset.f_clear==false) return false;

  varset.p_quab=new Dipole::Branch(ban);    //Copying!!
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);
  varset.p_root=root;

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

  varset.f_clear=false;

  m_memo=m_name;

  return true;

}





const bool Chain::StartInit(const Dipole::Branch&) {///////////////////////////
  assert(0);
  return false;
}
const bool Chain::StartInit(const Dipole::Glubranch&) {////////////////////////
  assert(0);
  return false;
}
const bool Chain::AddToInit(const Dipole::Glubranch&, const double&) {/////////
  assert(0);
  return false;
}
const bool Chain::FinishInit(const Dipole::Antibranch&, const double&) {///////
  assert(0);
  return false;
}
const bool Chain::FinishInit(const Dipole::Glubranch&, const double&) {////////
  assert(0);
  return false;
}



//=============================================================================



void Chain::Copy(const Chain& cha, const Type type) {

  //Later: Additional Chain initialization methods will cause further
  //modifications here!

  //Only two types are now possible: line and ring.
  //A line has one dipole at least.
  //A ring has two dipoles at least and the root gluon always resides at the
  //first place of the Glubranch list.

  varset.f_clear=false;
  varset.f_active=cha.varset.f_active;
  varset.m_k2tlast=cha.varset.m_k2tlast;

  if(type==line) {
    assert(!varset.p_quab); assert(!varset.p_atib);
    varset.p_quab=new Dipole::Branch(*cha.varset.p_quab);
    varset.p_atib=new Dipole::Antibranch(*cha.varset.p_atib);
    assert(varset.p_quab); assert(varset.p_atib);
    varset.p_1glu=cha.varset.p_1glu;    //In case a ring has been split.
  } else {
    assert(varset.l_glub.empty()); assert(!varset.p_1glu);
    varset.l_glub.push_front
      (new Dipole::Glubranch(*cha.varset.l_glub.front()));
    varset.p_1glu=varset.l_glub.front();
    assert(varset.p_1glu);
  }

  Dipole* dip=NULL;
  Dipole::Glubranch* preg=NULL, * botg=NULL;
  list<Dipole*>::const_iterator dit=cha.varset.l_dip.begin();
  list<Dipole*>::const_iterator dot=cha.varset.l_dip.end(); --dot;

  //First Dipole.
  if(dit!=dot) {
    Dipole_Particle* bo=(*dit)->GetBotBranchPointer().operator->();
    assert(bo->OrgType()==Nil);
    botg=static_cast<Dipole::Glubranch*>(bo);
    preg=botg=new Dipole::Glubranch(*botg); assert(botg);
    varset.l_glub.push_back(botg);
    if(type==line) dip=new Dipole(*varset.p_quab,*botg);
    else           dip=new Dipole(*varset.l_glub.front(),*botg);
    assert(dip);
    varset.l_dip.push_back(dip);
    if(*dit==cha.varset.p_root) varset.p_root=dip;
    ++dit;
  }
  //Inbetween Dipoles.
  for(; dit!=dot; ++dit) {
    Dipole_Particle* bo=(*dit)->GetBotBranchPointer().operator->();
    assert(bo->OrgType()==Nil);
    botg=static_cast<Dipole::Glubranch*>(bo);
    botg=new Dipole::Glubranch(*botg); assert(botg);
    varset.l_glub.push_back(botg);
    dip=new Dipole(*preg,*botg); assert(dip);
    varset.l_dip.push_back(dip);
    if(*dit==cha.varset.p_root) varset.p_root=dip;
    preg=botg;
  }
  //Last Dipole.
  if(!dip) {
    dip=new Dipole(*varset.p_quab,*varset.p_atib); assert(dip);
    varset.l_dip.push_back(dip);
    varset.p_root=dip;
  } else {
    if(type==line) dip=new Dipole(*preg,*varset.p_atib);
    else           dip=new Dipole(*preg,*varset.l_glub.front());
    assert(dip);
    varset.l_dip.push_back(dip);
    if(*dot==cha.varset.p_root) varset.p_root=dip;
  }

  if(cha.CheckMomentumConservation(varset.m_momentum)) {
    varset.m_mass=cha.varset.m_mass;
    varset.m_invmass=cha.varset.m_invmass;
  } else {
    varset.f_active=Blocked;
    if(varset.m_momentum[0]<0.0)
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Warning: Total energy is negative!\n"<<endl;
    varset.m_invmass=varset.m_momentum.Abs2();
    if(varset.m_invmass<0.0) {
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Warning: Negative invariant mass!\n"<<endl;
      varset.m_mass=-1*sqrt(-1*varset.m_invmass);
    }
    else varset.m_mass=sqrt(varset.m_invmass);
  }

}





const ATOOLS::Vec4D& Chain::UpdateMomentum(double k, const ATOOLS::Vec4D& p) {
  varset.f_clear=false;
  varset.m_momentum+=k*p;
  varset.m_invmass=varset.m_momentum.Abs2();
  if(varset.m_invmass<-1.0e-10/*0.0*/) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
	<<"Negative invariant mass ("<<varset.m_invmass<<") !\n"<<endl;
    //varset.f_active=Blocked;
    varset.m_mass=-1*sqrt(-1*varset.m_invmass);
  }
  else varset.m_mass=sqrt(varset.m_invmass);
  return varset.m_momentum;
}



//=============================================================================





//eof
