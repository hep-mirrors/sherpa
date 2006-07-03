//bof
//Version: 4 ADICIC++-0.0/2006/06/01

//Implementation of Dipole.H.


#ifdef __GNUC__
#if __GNUC__ >2
#include <ios>
#endif
#endif
#include <iomanip>
#include <string>
#include <sstream>
#include "Message.H"
#include "Dipole.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Dipole.exa.cc"





//=============================================================================



ostream& ADICIC::operator<<(ostream& ost, const ADICIC::Dipole& dip) {
  ost<<om::bold<<"Dipole "<<dip.m_name<<": "
     <<dip.p_top->Flav()<<" | "
     <<dip.p_bot->Flav()<<om::reset<<"    type=";
  switch(dip.m_type) {
  case -9999: ost<<"incorrect"; break;
  default   : ost<<dip.m_type;
  }
  if(dip.p_top->Incoming()) ost<<"    i"; else ost<<"    f";
  if(dip.p_bot->Incoming()) ost<<"i"; else ost<<"f";
  ost<<"    origs:"<<dip.p_top->OrgType()<<"|"<<dip.p_bot->OrgType();
  ost<<"    hdl:"<<bool(dip.p_hdl)<<"    ph:"<<dip.PointerHandling();
  ost<<"    state:";
  switch(dip.f_active) {
  case -1 : ost<<"blocked,"<<dip.m_nchg; break;
  case 0  : ost<<"off,"<<dip.m_nchg; break;
  default : ost<<"on,"<<dip.m_nchg;
  }
  ost<<"    spico:"<<dip.f_spico;
  ost<<endl<<setiosflags(ios::left)<<"  "
     <<"mass :"<<setw(12)<<dip.m_mass
     <<"sqrm :"<<setw(12)<<dip.m_invmass
     <<"facsc:"<<setw(12)<<dip.m_fc2
     <<"P="<<dip.m_momentum;
  string st1, st2, st3, st4; stringstream cv1, cv2, cv3, cv4;
  cv1<<dip.m_memory; cv2<<dip.m_copy; cv1>>st1; cv2>>st2;
  cv3<<dip.p_top->Name; cv4<<dip.p_bot->Name; cv3>>st3; cv4>>st4;
  st1="(memo:"+st1+")";
  st2="[copy:"+st2+"]";
  st3+="]"; st4+="]";
  ost<<endl<<"  pscal:"<<setw(12)<<dip.m_p2t
     <<"bscal:"<<setw(12)<<dip.m_k2t
     <<"escal:"<<setw(12)<<dip.m_l2t
     <<"p="<<dip.p_top->Momentum()<<"\t"<<dip.p_top->Momentum().Abs2();
  ost<<endl<<"  "<<setw(18)<<st1<<setw(18)<<st2
     <<"br["<<setw(6)<<st3<<"ab["<<setw(6)<<st4
     <<"q="<<dip.p_bot->Momentum()<<"\t"<<dip.p_bot->Momentum().Abs2();
  ost<<resetiosflags(ios::left);
  return ost;
}



//=============================================================================



void Dipole::Gate::operator()(Dipole* dipo,
			      Dipole_Particle* Od,
			      Dipole_Particle* Nw, bool b) const {

  if(dipo->p_top==Od) {
    dipo->p_top=Nw;
    if(b) dipo->f_top=true;
#ifdef DIPOLE_OUTPUT
    cout<<"    confirm TopBranch re-linking for Dipole "<<dipo->m_name
	<<"."<<endl;
#endif
    return;
  }
  if(dipo->p_bot==Od) {
    dipo->p_bot=Nw;
    if(b) dipo->f_bot=true;
#ifdef DIPOLE_OUTPUT
    cout<<"    confirm BotBranch re-linking for Dipole "<<dipo->m_name
	<<"."<<endl;
#endif
    return;
  }
  cerr<<"\nBug: Dipole_Particle does not belong to the Dipole!\n";
  assert(0);

}



//=============================================================================



int Dipole::s_count=0;    //So far not faked since there is no static Dipole.
int Dipole::s_maxcount=0;
const int& Dipole::InStore=Dipole::s_count;



//=============================================================================



Dipole::Dipole()
  : f_top(true), f_bot(true),
    m_name(++s_maxcount), m_copy(0), m_nchg(0),
    m_memory(0), f_active(Blocked),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect), f_spico(false),
    m_p2t(0.0), m_k2t(0.0), m_l2t(0.0), m_fc2(0.0),
    m_mass(0.0), m_invmass(0.0), m_momentum(Vec4D()),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  ++s_count;

  p_top=new Glubranch(); assert(p_top);
  p_bot=new Glubranch(); assert(p_bot);

  UpdateType();    //Increments the StateNumber by one.
  //UpdateMass();    //Is not necessary because nil Vec4Ds are used.

  AddDipoleToTowers();

}





Dipole::Dipole(const Dipole& dip, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), m_copy(dip.m_name), m_nchg(1),
    m_memory(dip.m_memory), f_active(dip.f_active),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(dip.m_type), f_spico(dip.f_spico),
    m_p2t(dip.m_p2t), m_k2t(dip.m_k2t), m_l2t(dip.m_l2t), m_fc2(dip.m_fc2),
    m_mass(dip.m_mass), m_invmass(dip.m_invmass), m_momentum(dip.m_momentum),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  //Due to "towering" type and momentum of Dipole dip are already up-to-date.

  ++s_count;

  if(f_top) {    //Flag f_bot will not give any new information.
    if(dip.p_top->OrgType()==Positive)
      p_top=new Branch(*static_cast<Branch*>(dip.p_top));
    else p_top=new Glubranch(*static_cast<Glubranch*>(dip.p_top));
    assert(p_top);
    if(dip.p_bot->OrgType()==Negative)
      p_bot=new Antibranch(*static_cast<Antibranch*>(dip.p_bot));
    else p_bot=new Glubranch(*static_cast<Glubranch*>(dip.p_bot));
    assert(p_bot);
  } else {
    if(dip.p_top->OrgType()==Positive)
      p_top=static_cast<Branch*>(dip.p_top);
    else p_top=static_cast<Glubranch*>(dip.p_top);
    if(dip.p_bot->OrgType()==Negative)
      p_bot=static_cast<Antibranch*>(dip.p_bot);
    else p_bot=static_cast<Glubranch*>(dip.p_bot);
  }

  //UpdateType(); UpdateMass();    //Is therefore not necessary.

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Branch& ban, Dipole::Antibranch& ati,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), m_copy(0), m_nchg(-1),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  ++s_count;

  if(f_top) {    //Flag f_bot will not give any new information.
    p_top=new Branch(ban);
    p_bot=new Antibranch(ati);
  } else {
    p_top=&ban;
    p_bot=&ati;
  }
  assert(p_top);
  assert(p_bot);

  UpdateType();
  UpdateMass();    //Increments the StateNumber by two altogether.

  f_spico=true; m_p2t=m_k2t=m_l2t=m_fc2=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Branch& ban, Dipole::Glubranch& glu,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), m_copy(0), m_nchg(-1),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  ++s_count;

  if(f_top) {    //Flag f_bot will not give any new information.
    p_top=new Branch(ban);
    p_bot=new Glubranch(glu);
  } else {
    p_top=&ban;
    p_bot=&glu;
  }
  assert(p_top);
  assert(p_bot);

  UpdateType();
  UpdateMass();

  f_spico=true; m_p2t=m_k2t=m_l2t=m_fc2=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Glubranch& glu, Dipole::Antibranch& ati,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), m_copy(0), m_nchg(-1),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  ++s_count;

  if(f_top) {    //Flag f_bot will not give any new information.
    p_top=new Glubranch(glu);
    p_bot=new Antibranch(ati);
  } else {
    p_top=&glu;
    p_bot=&ati;
  }
  assert(p_top);
  assert(p_bot);

  UpdateType();
  UpdateMass();

  f_spico=true; m_p2t=m_k2t=m_l2t=m_fc2=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Glubranch& glut, Dipole::Glubranch& glub,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), m_copy(0), m_nchg(-1),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect),
    Name(m_name), CopyOf(m_copy), StateNumber(m_nchg) {

  if(&glut==&glub && !f_top) {
    cerr<<"\nError: TopBranch and BotBranch point to a single gluon!\n";
    assert(0);
  }

  ++s_count;

  if(f_top) {    //Flag f_bot will not give any new information.
    p_top=new Glubranch(glut);
    p_bot=new Glubranch(glub);
  } else {
    p_top=&glut;
    p_bot=&glub;
  }
  assert(p_top);
  assert(p_bot);

  UpdateType();
  UpdateMass();

  f_spico=true; m_p2t=m_k2t=m_l2t=m_fc2=m_invmass;

  AddDipoleToTowers();

}





Dipole::~Dipole() {

  if(p_hdl) {    //Check for Dipole_Handler.
    cerr<<"\nMethod: ADICIC::Dipole::~Dipole(): "
	<<"Warning: Detaching Dipole_Handler from Dipole!\n"<<endl;
    if(p_hdl->IsDockedAt(*this)==false) {
      cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
      assert(p_hdl->IsDockedAt(*this));
    }
    Dipole_Handler* phand=p_hdl;
    p_hdl=NULL;
    phand->DetachDipole(this);
  }

  RemoveDipoleFromTowers();    //To be sure.

#ifdef DIPOLE_OUTPUT
  cout<<"    destructing Dipole ["<<m_name<<"]"<<endl;
#endif

  if(f_top) {
    list<Dipole*>& toptow=AccessTopBranch();
    if(toptow.empty()) { if(p_top) delete p_top;}
    else {
      Test(toptow,this->p_top,"Top");
      Dipole* payee=toptow.front();
      if(payee->p_top==this->p_top) payee->f_top=true;
      else payee->f_bot=true;
    }
  }

  if(f_bot) {
    list<Dipole*>& bottow=AccessBotBranch();
    if(bottow.empty()) { if(p_bot) delete p_bot;}
    else {
      Test(bottow,this->p_bot,"Bot");
      Dipole* payee=bottow.front();
      if(payee->p_top==this->p_bot) payee->f_top=true;
      else payee->f_bot=true;
    }
  }

  p_top=NULL; p_bot=NULL;

  --s_count;
  if(m_name && m_name==s_maxcount) --s_maxcount;
  if(!s_count) s_maxcount=0;

}





Dipole& Dipole::operator=(const Dipole& dip) {

  if(this==&dip) return *this;

  m_copy=dip.m_name;
  ++m_nchg;
  m_memory=dip.m_memory;
  f_active=dip.f_active;

  m_type=dip.m_type;
  f_spico=dip.f_spico; 
  m_p2t=dip.m_p2t;
  m_k2t=dip.m_k2t;
  m_l2t=dip.m_l2t;
  m_fc2=dip.m_fc2;
  m_mass=dip.m_mass;
  m_invmass=dip.m_invmass;
  m_momentum=dip.m_momentum;

  this->RemoveDipoleFromTowers();

  if(f_top) {
    list<Dipole*>& totow=AccessTopBranch();
    if(totow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
    else {
      Test(totow,this->p_top,"Top");
      Dipole* payee=totow.front();
      if(payee->p_top==this->p_top) payee->f_top=true; else payee->f_bot=true;
    }
    if(dip.p_top->OrgType()==Positive)
      p_top=new Branch(*static_cast<Branch*>(dip.p_top));
    else p_top=new Glubranch(*static_cast<Glubranch*>(dip.p_top));
    assert(p_top);
  }
  else {
    if(dip.p_top->OrgType()==Positive) p_top=static_cast<Branch*>(dip.p_top);
    else p_top=static_cast<Glubranch*>(dip.p_top);
  }

  if(f_bot) {
    list<Dipole*>& botow=AccessBotBranch();
    if(botow.empty()) { if(p_bot) { delete p_bot; p_bot=NULL;}}
    else {
      Test(botow,this->p_bot,"Bot");
      Dipole* payee=botow.front();
      if(payee->p_top==this->p_bot) payee->f_top=true; else payee->f_bot=true;
    }
    if(dip.p_bot->OrgType()==Negative)
      p_bot=new Antibranch(*static_cast<Antibranch*>(dip.p_bot));
    else p_bot=new Glubranch(*static_cast<Glubranch*>(dip.p_bot));
    assert(p_bot);
  }
  else {
    if(dip.p_bot->OrgType()==Negative)
      p_bot=static_cast<Antibranch*>(dip.p_bot);
    else p_bot=static_cast<Glubranch*>(dip.p_bot);
  }

  AddDipoleToTowers();

  return *this;

}





const bool Dipole::operator==(const Dipole& dip) const {
  cerr<<"\nSorry :o( Method has not been implemented yet.\n";
  assert(0); exit(1);
}





void Dipole::Print() const {
  cout<<endl<<string(80,'=')<<endl;
  cout<<(*this)<<endl;
  cout<<"TopBranch: "<<p_top->Name;
  p_top->ShowParticle();
  cout<<"BotBranch: "<<p_bot->Name;
  p_bot->ShowParticle();
  cout<<string(80,'=')<<endl;
}





const Dipole::Type Dipole::UpdateType() {
  if(p_top->OrgType()==-1 || p_bot->OrgType()==1) {    //Cast protector!
    cerr<<"\nError: Demand will produce invalid Dipole type!\n";
    assert( (p_top->OrgType()==-1 || p_bot->OrgType()==1) == false );
  }
  ++m_nchg;
  short t = 10*p_top->Tag() + p_bot->Tag();
  //Does that sufficiently work? No.
  //Only if invalid Dipoles are really excluded.
#ifndef DIPOLE_STRICT_VERSION
  if(t==9 || t==10 || t==-1 || t==0); else {
    cerr<<"\nError: Algorithm used produces invalid Dipole type!\n";
    assert(t==9 || t==10 || t==-1 || t==0);
  }
  //#error
#endif
  if(t<0) t-=(1000*p_top->Incoming() + 100*p_bot->Incoming());
  else    t+=(1000*p_top->Incoming() + 100*p_bot->Incoming());
  m_type=Type(t);
  return m_type;
}





const double Dipole::UpdateMass() {
  ++m_nchg;
  if(p_top->Incoming()) m_momentum=(-1)*p_top->Momentum();
  else m_momentum=p_top->Momentum();
  if(p_bot->Incoming()) m_momentum+=(-1)*p_bot->Momentum();
  else m_momentum+=p_bot->Momentum();
  //switch(2*p_top->Incoming()+p_bot->Incoming()) {
  //case  2: m_momentum=p_top->Momentum()-p_bot->Momentum(); break;
  //case  1: m_momentum=p_bot->Momentum()-p_top->Momentum(); break;
  //default: m_momentum=p_top->Momentum()+p_bot->Momentum();
  //}
  m_invmass=m_momentum.Abs2();
  if(m_invmass<0.0) {
    //cerr<<"\nMethod: const double ADICIC::Dipole::UpdateMass(): "
    //    <<"Warning: Negative invariant mass!\n"<<endl;
    m_mass=-1*sqrt(-1*m_invmass);
    //The minus sign functions only as a flag.
  }
  else m_mass=sqrt(m_invmass);
  return m_mass;
}





void Dipole::RenewBranch(Branch& ban) {

  //This method is only concerned with p_top.

  if(p_top==&ban) return;

  list<Dipole*>& oldtow=AccessTopBranch();
  oldtow.remove(this);

  if(f_top) {
    if(oldtow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
    else {
      Test(oldtow,this->p_top,"Top");
      Dipole* payee=oldtow.front();
      if(payee->p_top==this->p_top) payee->f_top=true; else payee->f_bot=true;
    }
    p_top=new Branch(ban); assert(p_top);
  }
  else p_top=&ban;

  --m_nchg;
  UpdateType();    //That concerns only this dipole.
  UpdateMass();

  list<Dipole*>& newtow=AccessTopBranch();
  newtow.remove(this);
  newtow.push_front(this);

}





void Dipole::RenewBranch(Antibranch& ati) {

  //This method is only concerned with p_bot.

  if(p_bot==&ati) return;

  list<Dipole*>& oldtow=AccessBotBranch();
  oldtow.remove(this);

  if(f_bot) {
    if(oldtow.empty()) { if(p_bot) { delete p_bot; p_bot=NULL;}}
    else {
      Test(oldtow,this->p_bot,"Bot");
      Dipole* payee=oldtow.front();
      if(payee->p_top==this->p_bot) payee->f_top=true; else payee->f_bot=true;
    }
    p_bot=new Antibranch(ati); assert(p_bot);
  }
  else p_bot=&ati;

  --m_nchg;
  UpdateType();    //That concerns only this dipole.
  UpdateMass();

  list<Dipole*>& newtow=AccessBotBranch();
  newtow.remove(this);
  newtow.push_front(this);

}





void Dipole::RenewBranch(bool top, Glubranch& glu) {

  if(top) {    //Then method is concerned with p_top.

    if(p_top==&glu) return;

    list<Dipole*>& oldtow=AccessTopBranch();
    oldtow.remove(this);

    if(f_top) {
      if(oldtow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
      else {
	Test(oldtow,this->p_top,"Top");
	Dipole* payee=oldtow.front();
	if(payee->p_top==this->p_top) payee->f_top=true;
	else payee->f_bot=true;
      }
      p_top=new Glubranch(glu); assert(p_top);
    }
    else p_top=&glu;

    list<Dipole*>& newtow=AccessTopBranch();
    newtow.remove(this);
    newtow.push_front(this);

  } else {    //Then method is concerned with p_bot.

    if(p_bot==&glu) return;

    list<Dipole*>& oldtow=AccessBotBranch();
    oldtow.remove(this);

    if(f_bot) {
      if(oldtow.empty()) { if(p_bot) { delete p_bot; p_bot=NULL;}}
      else {
	Test(oldtow,this->p_bot,"Bot");
	Dipole* payee=oldtow.front();
	if(payee->p_top==this->p_bot) payee->f_top=true;
	else payee->f_bot=true;
      }
      p_bot=new Glubranch(glu); assert(p_bot);
    }
    else p_bot=&glu;

    list<Dipole*>& newtow=AccessBotBranch();
    newtow.remove(this);
    newtow.push_front(this);

  }

  --m_nchg;
  UpdateType();    //That concerns only this dipole.
  UpdateMass();

}





const bool Dipole::RenewBranch(bool top, const Dipole& dip, bool diptop) {
  cerr<<"\nSorry :o( Method has not been implemented yet. ";
  cerr<<"Not sure if it is really needed.\n";
  assert(0); exit(1);
  return false;
  //if(this!=&dip)    //Maybe dropping that to allow for cross assignments.
  //gg check
}



//=============================================================================



void Dipole::Test(const std::list<Dipole*>& L, const Dipole_Particle* P,
		  const std::string& s) const {
  for(std::list<Dipole*>::const_iterator it=L.begin(); it!=L.end(); ++it) {
    if( (*it)->p_top==P ) {
      if( (*it)->f_top ) {
	cerr<<"\nBug: "<<s
	    <<"Branch physically belongs to several Dipoles!\n";
	assert((*it)->f_top==false);}
    } else {
      if( (*it)->p_bot==P ) {
	if( (*it)->f_bot ) {
	  cerr<<"\nBug: "<<s
	      <<"Branch physically belongs to several Dipoles!\n";
	  assert((*it)->f_bot==false);}
      } else {
	cerr<<"\nBug: Dipole_Particle's tower carries external Dipole!";
	cerr<<endl; assert(0);
      }
    }
  }
}



//=============================================================================





//eof
