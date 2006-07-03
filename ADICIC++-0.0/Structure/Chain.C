//bof
//Version: 4 ADICIC++-0.0/2006/06/14

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
     <<cha.Source<<" =="<<om::reset<<"\n";
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
  ost<<"    Chain handler: "<<cha.IsHandled();
  ost<<"    #(p|d|cp): "<<cha.ParticleNumber()<<"|"
     <<cha.DipoleNumber()<<"|"<<cha.CorrParticleNumber()<<"\n";
  ost<<setiosflags(ios::left)
     <<"EvoScale="<<setw(13)<<cha.LastScale()
     <<"Fa"<<dpa.evo.FactScaleType()<<"Scale="<<setw(13)<<cha.FactScale()
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
      ost<<om::bold<<"Glubranch(es):"<<om::reset<<"\n";
    for(list<Dipole::Glubranch*>::const_iterator
	  gut=cha.GlubranchPointerList().begin();
	gut!=cha.GlubranchPointerList().end(); ++gut) {
      if(*gut==cha.FirstGlubranchPointer())
	ost<<(*gut)->Parton<<"  "<<om::red<<"g_root"<<om::reset<<"\n";
      else
	ost<<(*gut)->Parton<<"\n";
    }
    if(!cha.CorrParticlePointerList().empty())
      ost<<om::bold<<"Correlated Particle(s):"<<om::reset<<"\n";
    for(list<Particle*>::const_iterator
	  wat=cha.CorrParticlePointerList().begin();
	wat!=cha.CorrParticlePointerList().end(); ++wat) ost<<(**wat)<<"\n";
    ost<<om::bold<<" "<<string(118,'-')<<om::reset<<"\n";
    ost<<"Number of  branches = "<<cha.ParticleNumber()<<"\n";
    ost<<"Number of   dipoles = "<<cha.DipoleNumber()<<"\n";
    if(!cha.CorrParticlePointerList().empty())
      ost<<"Number of corr ptcs = "<<cha.CorrParticleNumber()<<"\n";
  }
  ost<<om::bold<<string(114,'=')<<" "<<cha.Name<<" | "
     <<cha.Source<<" =="<<om::reset;
  return ost;
}



//=============================================================================



const Chain::Initiator::EpEm Chain::Initiator::epem={};



//=============================================================================



Chain::Kernel::Kernel() { Init();}



const bool Chain::Kernel::IsInit() const {
  return (f_active==Blocked && m_k2tlast==0.0 && m_mufac.empty() &&
	  p_root==NULL && p_hdl==NULL &&
	  p_quab==NULL && p_atib==NULL && p_1glu==NULL &&
	  m_momentum[0]==0.0 && m_momentum[1]==0.0 &&
	  m_momentum[2]==0.0 && m_momentum[3]==0.0 &&
	  l_glub.empty() && l_dip.empty() && l_corr.empty());
}



void Chain::Kernel::Init() {
  f_active=Blocked;
  m_k2tlast=m_mass=m_invmass=0.0; m_momentum=Vec4D();
  m_mufac.clear();
  p_root=NULL;
  p_hdl=NULL;
  p_quab=NULL; p_atib=NULL; p_1glu=NULL;
  l_glub.clear(); l_dip.clear(); l_corr.clear();
}



void Chain::Kernel::Destruct() {
  m_mufac.clear();
  for(list<Particle*>::iterator pat=l_corr.begin(); pat!=l_corr.end(); ++pat)
    if(*pat) delete (*pat);
  for(list<Dipole*>::iterator dit=l_dip.begin(); dit!=l_dip.end(); ++dit)
    if(*dit) delete (*dit);
  for(list<Dipole::Glubranch*>::iterator git=l_glub.begin();
      git!=l_glub.end(); ++git) if(*git) delete (*git);
  if(p_quab) delete p_quab;
  if(p_atib) delete p_atib;
  l_corr.clear();
  l_dip.clear();
  l_glub.clear();
}



//=============================================================================



bool Chain::s_print=false;
bool& Chain::MoreCout=Chain::s_print;

size_t Chain::s_count=0;
size_t Chain::s_maxcount=0;
const size_t& Chain::InStore=Chain::s_count;
const size_t& Chain::MaxCount=Chain::s_maxcount;



//=============================================================================



Chain::Chain()
  : m_name(++s_maxcount), m_memo(0), varset(), l_fill(),
    Name(m_name), Source(m_memo) {
  ++s_count;
  //This is a clear chain.
}





Chain::Chain(const Chain& cha)
  : m_name(++s_maxcount), m_memo(0), varset(), l_fill(),
    Name(m_name), Source(m_memo) {
  ++s_count;
  //Later: Additional Chain initialization methods will cause further
  //modifications here!
  Type type=cha.ChainType();
  if(type==incorrect) return;//////////////////////////////////////////////////
  //Chain cha is not empty here!
  //Do not copy cha.l_fill!
  this->Copy(cha,type);
  m_memo=cha.m_memo;
}





Chain::Chain(const Dipole::Branch& ban, const Dipole::Antibranch& ati,
	     const Initiator::EpEm)
  : m_name(++s_maxcount), m_memo(0), varset(), l_fill(),
    Name(m_name), Source(m_memo) {

  //Only EpEm -> FF chain!
  assert(ban.Incoming()==false && ati.Incoming()==false);

  ++s_count;

  varset.p_quab=new Dipole::Branch(ban);    //Copying!!
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);
  varset.p_root=root;

  varset.m_k2tlast=root->ProdScale();
  varset.m_mufac.resize(fascat::stop,varset.m_k2tlast);
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

}





Chain::Chain(const Dipole::Glubranch& glut, const Dipole::Glubranch& glub,
	     const Initiator::EpEm)
  : m_name(++s_maxcount), m_memo(0), varset(), l_fill(),
    Name(m_name), Source(m_memo) {

  //Only EpEm -> FF chain!
  assert(glut.Incoming()==false && glub.Incoming()==false);

  ++s_count;

  if(&glut==&glub) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Gluon Chain initiation from a single gluon!\n"<<endl;
    //exit(1);
  }

  Dipole::Glubranch* top=new Dipole::Glubranch(glut);    //Copying!!
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
  varset.m_mufac.resize(fascat::stop,varset.m_k2tlast);
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

  for(list<Particle*>::iterator pat=l_fill.begin(); pat!=l_fill.end(); ++pat)
    if(*pat) delete (*pat);

  varset.Destruct();

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
  //The message will come from Clear() itself.
  //Do not copy cha.l_fill, i.e. l_fill remains as it is.
  this->Copy(cha,type);
  m_memo=cha.m_memo;
  return *this;
}





//const bool Chain::operator==(const Chain& cha) const {
//  cerr<<"\nSorry :o( Method has not been implemented.\n";
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





const unsigned Chain::INumber() const {
  if(this->IsEmpty()) return 0;
  unsigned in=0;
  if(varset.p_quab) if(varset.p_quab->Incoming()) ++in;
  if(varset.p_atib) if(varset.p_atib->Incoming()) ++in;
  for(list<Dipole::Glubranch*>::const_iterator gut=varset.l_glub.begin();
      gut!=varset.l_glub.end(); ++gut) {
    assert(*gut);
    if((*gut)->Incoming()) ++in;
  }
  assert(in<3);
  if(in==2) assert(!varset.l_corr.empty());
  return in;
}





const bool Chain::CheckMomentumConservation(Vec4D& sum) const {
  static size_t viol=0;
  sum=Vec4D();
  unsigned in=0;
  if(varset.p_quab) {
    if(varset.p_quab->Incoming()) { sum-=varset.p_quab->Momentum(); ++in;}
    else sum+=varset.p_quab->Momentum();
  }
  if(varset.p_atib) {
    if(varset.p_atib->Incoming()) { sum-=varset.p_atib->Momentum(); ++in;}
    else sum+=varset.p_atib->Momentum();
  }
  for(list<Dipole::Glubranch*>::const_iterator gut=varset.l_glub.begin();
      gut!=varset.l_glub.end(); ++gut) {
    assert(*gut);
    if((*gut)->Incoming()) { sum-=(*gut)->Momentum(); ++in;}
    else sum+=(*gut)->Momentum();
  }
  assert(in<3);
  Vec4D test=Vec4D();
  if(in==2) {
    //Correlated particles are potentially changed only in chains
    //containing an II dipole.
    for(list<Particle*>::const_iterator pat=varset.l_corr.begin();
	pat!=varset.l_corr.end(); ++pat) {
      assert(*pat); test+=(*pat)->Momentum();
    }
    test+=sum;
    for(char icorr=0; icorr<4; ++icorr) { if(dabs(test[icorr])>1.0e-9)
      cout<<__PRETTY_FUNCTION__<<": vec-conservation violation="
	  <<(++viol)<<endl;}    //assert(dabs(test[icorr])<1.0e-8);
  }
  if(sum==varset.m_momentum) return true;
  test=sum-varset.m_momentum;
  for(char i=0; i<4; ++i) if(dabs(test[i])>1.0e-10) return false;
  return true;
}





const boolint Chain::ExtractPartons(Particle_List& parlist) const {

  //The chain entries are appended to the parlist.

  Type type=this->ChainType();
  if(type==incorrect) return boolint();    //(false,0)!

  list<Dipole*>::const_iterator dit=varset.l_dip.begin();
  Particle* par=new Particle( (*dit)->GetTopBranchPointer()->Parton );
  if((*dit)->GetTopBranchPointer()->Incoming()) par->SetInfo('i');
  else par->SetInfo('f');
  parlist.push_back(par);
  int minnum=par->Number();
  Particle_List::iterator pit=parlist.end(); --pit;
  ++dit;
  for(; dit!=varset.l_dip.end(); ++dit) {
    par=new Particle( (*dit)->GetTopBranchPointer()->Parton );
    if((*dit)->GetTopBranchPointer()->Incoming()) par->SetInfo('i');
    else par->SetInfo('f');
    parlist.push_back(par);
    if(par->Number()<minnum) minnum=par->Number();
  }

  if(type==line) {
    --dit;
    par=new Particle( (*dit)->GetBotBranchPointer()->Parton );
    if((*dit)->GetBotBranchPointer()->Incoming()) par->SetInfo('i');
    else par->SetInfo('f');
    parlist.push_back(par);
    if(par->Number()<minnum) minnum=par->Number();

    const size_t num=parlist.size(); assert(num>1);
    if(parlist.size()>2) {
      switch(this->INumber()%2) {
      case  0: {
	bool z=(*pit)->Info()=='i';
	unsigned flow=(*pit)->GetFlow(1+z); assert((*pit)->GetFlow(2-z)==0);
	if(flow==0) { (*pit)->SetFlow(1+z,-1); flow=(*pit)->GetFlow(1+z);}
	++pit;
	Particle_List::iterator pot=parlist.end();
	--pot;
	for(; pit!=pot; ++pit) {
	  z=(*pit)->Info()=='i';
	  (*pit)->SetFlow(2-z,flow);
	  (*pit)->SetFlow(1+z,-1);
	  flow=(*pit)->GetFlow(1+z);
	}
	z=(*pit)->Info()=='i';
	(*pit)->SetFlow(2-z,flow); assert((*pot)->GetFlow(1+z)==0);
      } break;
      case  2: {
	unsigned doflow=(*pit)->GetFlow(2); assert((*pit)->GetFlow(1)==0);
	if(doflow==0) { (*pit)->SetFlow(2,-1); doflow=(*pit)->GetFlow(2);}
	++pit;
	Particle_List::iterator pot=parlist.end();
	--pot;
	for(; pit!=pot; ++pit) {
	  (*pit)->SetFlow(2,doflow);
	  (*pit)->SetFlow(1,-1);
	  doflow=(*pit)->GetFlow(1);
	}
	(*pit)->SetFlow(1,doflow); assert((*pot)->GetFlow(2)==0);
      } break;
      default: assert(0);
      }
    } else {
      switch(this->INumber()) {
      case  0: {
	unsigned upflow=(*pit)->GetFlow(1); assert((*pit)->GetFlow(2)==0);
	if(upflow==0) { (*pit)->SetFlow(1,-1); upflow=(*pit)->GetFlow(1);}
	++pit; assert((*pit)->GetFlow(1)==0);
	if((*pit)->GetFlow(2)!=upflow) (*pit)->SetFlow(2,upflow);
      } break;
      case  2: {
	unsigned doflow=(*pit)->GetFlow(2); assert((*pit)->GetFlow(1)==0);
	if(doflow==0) { (*pit)->SetFlow(2,-1); doflow=(*pit)->GetFlow(2);}
	++pit; assert((*pit)->GetFlow(2)==0);
	if((*pit)->GetFlow(1)!=doflow) (*pit)->SetFlow(1,doflow);
      } break;
      default: assert(0);
      }
    }
  }
  else {
    assert(0);
    //type==ring ... No flow setting so far.
  }

  if(varset.l_corr.empty()); else {
    list<Particle*>::const_iterator cit=varset.l_corr.begin();
    for(; cit!=varset.l_corr.end(); ++cit) {
      par=new Particle(**cit);
      parlist.push_back(par);
    }
  }

#ifdef CHAIN_OUTPUT
  cout<<parlist<<endl;
#endif

  return boolint(true,minnum);

}





const bool Chain::Clear() {
  if(varset.p_hdl) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Clearing a Chain "
	<<"being manipulated by a Chain_Handler is not permitted!\n"<<endl;
    return false;
  }
  //Leave l_fill untouched.
  varset.Destruct();
  varset.Init();    //This is a clear chain right now.
  m_memo=0;
  return true;
}





const bool Chain::AddCorrParticleToInit(const Particle& par) {
  Flavour flav=par.RefFlav();
  if(flav.IsLepton() ||
     flav.Strong() ||
     flav.IsPhoton());
  else return false;
  l_fill.push_back(new Particle(par));    //Copying!!
  assert(l_fill.back());
  l_fill.back()->SetInfo('c');
  return true;
}





const bool Chain::Initialize(const Dipole::Branch& ban,
			     const Dipole::Antibranch& ati,
			     double scale) {

  if(!IsClear()) return false;

  //So we have m_memo==0.

  varset.p_quab=new Dipole::Branch(ban);    //Copying!!
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);
  varset.p_root=root;

  varset.m_mufac.resize(fascat::stop,root->ProdScale());

  if(scale) root->SetEvolScales(scale);

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();

  Vec4D test=root->TotP();
  list<Particle*>::iterator pat=l_fill.begin();
  if(root->IsFF()) {
    if(varset.m_invmass<.0 || varset.m_momentum[0]<.0) varset.f_active=Blocked;
    //Assume corrs = f+f = dip (momenta in in=out definition (flow)).
    if(!l_fill.empty()) {
      for(; pat!=l_fill.end(); ++pat) test-=(*pat)->Momentum();
      for(char ffco=0; ffco<4; ++ffco) { if(dabs(test[ffco])>1.0e-10) {
	cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Proposal for "
	    <<"correlated particles is ignored due to momentum mismatch!\n\n";
	return true;
      }}
      varset.l_corr.splice(varset.l_corr.begin(),l_fill);
      assert(l_fill.empty());
    }
    return true;
  } else {
    if(root->IsII()) {
      if(varset.m_invmass<0.0) varset.f_active=Blocked;
      //Assume -i-i + corrs = 0 = dip + corrs.
      if(!l_fill.empty()) {
	bool flag=true;
	for(; pat!=l_fill.end(); ++pat) test+=(*pat)->Momentum();
	for(char iico=0; iico<4; ++iico) {
	  if(dabs(test[iico])>1.0e-10) flag=false; break;}
	if(flag) {
	  varset.l_corr.splice(varset.l_corr.begin(),l_fill);
	  assert(l_fill.empty());
	  return true;
	}
      }
      //Creating the residual particle!!
      varset.l_corr.push_back(new Particle(-7,Flavour(kf::cluster),
					   (-1.0)*root->TotP()));
      assert(varset.l_corr.back());
      varset.l_corr.back()->SetInfo('r');
      return true;
    } else {
      if(varset.m_invmass>0.0) varset.f_active=Blocked;
      //Assume dip = f-i = corri - corrs and 1st particle of l_fill is corri.
      if(!l_fill.empty()) {
	test-=(*pat)->Momentum(); ++pat;
	for(; pat!=l_fill.end(); ++pat) test+=(*pat)->Momentum();
	for(char fico=0; fico<4; ++fico) { if(dabs(test[fico])>1.0e-10) {
	  cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Proposal for "
	      <<"correlated particles is ignored! Momentum mismatch!\n\n";
	  return true;
	}}
	varset.l_corr.splice(varset.l_corr.begin(),l_fill);
	assert(l_fill.empty());
      }
      return true;
    }
  }

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

  varset.f_active=cha.varset.f_active;
  varset.m_k2tlast=cha.varset.m_k2tlast;
  varset.m_mufac=cha.varset.m_mufac;    //Real copies are created.

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
    dip->ApplySettingsOf(**dit);
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
    dip->ApplySettingsOf(**dit);
    if(*dit==cha.varset.p_root) varset.p_root=dip;
    preg=botg;
  }
  //Last Dipole.
  if(!dip) {
    dip=new Dipole(*varset.p_quab,*varset.p_atib); assert(dip);
    varset.l_dip.push_back(dip);
    dip->ApplySettingsOf(**dit);
    varset.p_root=dip;
  } else {
    if(type==line) dip=new Dipole(*preg,*varset.p_atib);
    else           dip=new Dipole(*preg,*varset.l_glub.front());
    assert(dip);
    varset.l_dip.push_back(dip);
    dip->ApplySettingsOf(**dit);
    if(*dot==cha.varset.p_root) varset.p_root=dip;
  }

  //Copy the part of correlated particles.
  varset.l_corr=cha.varset.l_corr;
  for(list<Particle*>::iterator pit=varset.l_corr.begin();
      pit!=varset.l_corr.end(); ++pit) {
    *pit=new Particle(**pit); assert(*pit);
  }

  //Test 4-momenta.
  if(cha.CheckMomentumConservation(varset.m_momentum)) {
    varset.m_mass=cha.varset.m_mass;
    varset.m_invmass=cha.varset.m_invmass;
  } else {
    //Probably, this part is not really needed.
    varset.f_active=Blocked;
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Momentum conservation check failed!\n"<<endl;
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

  //Translate source information.
  //Relies on: preadjusted source is source of the dipole copied from.
  dit=cha.varset.l_dip.begin();
  dot=varset.l_dip.begin();
  for(; dit!=cha.varset.l_dip.end(); ++dit) {
    int oname=(*dit)->Name;
    int nname=(*dot)->Name;
    for(list<Dipole*>::iterator ot=varset.l_dip.begin();
	ot!=varset.l_dip.end(); ++ot)
      if((*ot)->Source()) if((*ot)->Source()==oname) (*ot)->SetSource()=nname;
    ++dot;
  }

}





const Vec4D& Chain::UpdateMomentum(double k, const Vec4D& p) {
  varset.m_momentum+=k*p;
  varset.m_invmass=varset.m_momentum.Abs2();
  if(varset.m_invmass<0.0) {
    //cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
    //    <<"Negative invariant mass ("<<varset.m_invmass<<") !\n"<<endl;
    //varset.f_active=Blocked;//???????????????????????????????????????????????
    varset.m_mass=-1*sqrt(-1*varset.m_invmass);
    //The minus sign functions only as a flag.
  }
  else varset.m_mass=sqrt(varset.m_invmass);
  return varset.m_momentum;
}



//=============================================================================





//eof
