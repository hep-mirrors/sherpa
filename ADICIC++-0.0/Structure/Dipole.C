//bof
//Version: 1 ADICIC++-0.0/2004/03/12

//Implementation of Dipole.H.



#include <ios>
#include <iomanip>
#include <string>
#include <sstream>
#include "Dipole.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Dipole.exa.cc"





//=============================================================================



ostream& ADICIC::operator<<(ostream& ost, const ADICIC::Dipole& dip) {
  ost<<"\e[1mDipole "<<dip.m_name<<": "
     <<dip.p_top->Flav()<<" | "
     <<dip.p_bot->Flav()<<"\e[0m   type:";
  switch(dip.m_type) {
  case -99: ost<<"incorrect"; break;
  default : ost<<dip.m_type;
  }
  ost<<"   origs:"<<dip.p_top->OrgType()<<"|"<<dip.p_bot->OrgType();
  ost<<"   hdl:"<<bool(dip.p_hdl);
  ost<<endl<<setiosflags(ios::left)<<"  "
     <<"mass :"<<setw(10)<<dip.m_mass
     <<"sqrm :"<<setw(12)<<dip.m_invmass
     <<"ph "<<setw(8)<<dip.PointerHandling()
     <<"P="<<dip.m_momentum;
  string st1, st2, st3, st4; stringstream cv1, cv2, cv3, cv4;
  cv1<<dip.m_memory; cv2<<dip.m_copy; cv1>>st1; cv2>>st2;
  cv3<<dip.p_top->Name; cv4<<dip.p_bot->Name; cv3>>st3; cv4>>st4;
  st1=string("(memo:")+st1+string(")");
  st2=string("[copy:")+st2+string("]");
  st3+=string("]"); st4+=string("]");
  ost<<endl<<"  state:"<<setw(10);
  switch(dip.f_active) {
  case -1 : ost<<"blocked"; break;
  case 0  : ost<<"off"; break;
  default : ost<<"on";
  }
  ost<<"scale:"<<setw(12)<<dip.m_k2t
     <<"br["<<setw(8)<<st3<<"p="<<dip.p_top->Momentum()
     <<"\t"<<dip.p_top->Momentum().Abs2();
  ost<<endl<<"  "<<setw(16)<<st1<<setw(18)<<st2
     <<"ab["<<setw(8)<<st4<<"q="<<dip.p_bot->Momentum()
     <<"\t"<<dip.p_bot->Momentum().Abs2();
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
    m_name(++s_maxcount), Name(m_name),
    m_copy(0), CopyOf(m_copy),
    m_memory(0), f_active(Blocked),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect), m_k2t(0.0),
    m_mass(0.0), m_invmass(0.0), m_momentum(Vec4D()) {

  ++s_count;

  p_top=new Glubranch(); assert(p_top);
  p_bot=new Glubranch(); assert(p_bot);

  UpdateType();
  //UpdateMass();    //Is not necessary because nil Vec4Ds are used.

  AddDipoleToTowers();

}





Dipole::Dipole(const Dipole& dip, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), Name(m_name),
    m_copy(dip.m_name), CopyOf(m_copy),
    m_memory(dip.m_memory), f_active(dip.f_active),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(dip.m_type), m_k2t(dip.m_k2t),
    m_mass(dip.m_mass), m_invmass(dip.m_invmass), m_momentum(dip.m_momentum) {

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
    m_name(++s_maxcount), Name(m_name),
    m_copy(0), CopyOf(m_copy),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect) {

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
  UpdateMass();

  m_k2t=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Branch& ban, Dipole::Glubranch& glu,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), Name(m_name),
    m_copy(0), CopyOf(m_copy),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect) {

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

  m_k2t=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Glubranch& glu, Dipole::Antibranch& ati,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), Name(m_name),
    m_copy(0), CopyOf(m_copy),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect) {

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

  m_k2t=m_invmass;

  AddDipoleToTowers();

}





Dipole::Dipole(Dipole::Glubranch& glut, Dipole::Glubranch& glub,
	       int source, bool phdl)
  : f_top(phdl), f_bot(phdl),
    m_name(++s_maxcount), Name(m_name),
    m_copy(0), CopyOf(m_copy),
    m_memory(source), f_active(On),
    p_top(NULL), p_bot(NULL), p_hdl(NULL),
    m_type(incorrect) {

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

  m_k2t=m_invmass;

  AddDipoleToTowers();

}





Dipole::~Dipole() {

  //Check for Dipole_Handler.
  if(p_hdl) {
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
      Test(toptow,this->p_top,string("Top"));
      Dipole* payee=toptow.front();
      if(payee->p_top==this->p_top) payee->f_top=true;
      else payee->f_bot=true;
    }
  }

  if(f_bot) {
    list<Dipole*>& bottow=AccessBotBranch();
    if(bottow.empty()) { if(p_bot) delete p_bot;}
    else {
      Test(bottow,this->p_bot,string("Bot"));
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
  m_memory=dip.m_memory;
  f_active=dip.f_active;

  m_type=dip.m_type;
  m_k2t=dip.m_k2t;
  m_mass=dip.m_mass;
  m_invmass=dip.m_invmass;
  m_momentum=dip.m_momentum;

  this->RemoveDipoleFromTowers();

  if(f_top) {
    list<Dipole*>& totow=AccessTopBranch();
    if(totow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
    else {
      Test(totow,this->p_top,string("Top"));
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
      Test(botow,this->p_bot,string("Bot"));
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





void Dipole::RenewBranch(Branch& ban) {

  //This method is only concerned with p_top.

  if(p_top==&ban && !f_top) return;

  list<Dipole*>& oldtow=AccessTopBranch();
  oldtow.remove(this);

  if(f_top) {
    if(oldtow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
    else {
      Test(oldtow,this->p_top,string("Top"));
      Dipole* payee=oldtow.front();
      if(payee->p_top==this->p_top) payee->f_top=true; else payee->f_bot=true;
    }
    p_top=new Branch(ban); assert(p_top);
  }
  else p_top=&ban;

  UpdateType();    //That concerns only this dipole.
  UpdateMass();

  list<Dipole*>& newtow=AccessTopBranch();
  newtow.remove(this);
  newtow.push_front(this);

}





void Dipole::RenewBranch(Antibranch& ati) {

  //This method is only concerned with p_bot.

  if(p_bot==&ati && !f_bot) return;

  list<Dipole*>& oldtow=AccessBotBranch();
  oldtow.remove(this);

  if(f_bot) {
    if(oldtow.empty()) { if(p_bot) { delete p_bot; p_bot=NULL;}}
    else {
      Test(oldtow,this->p_bot,string("Bot"));
      Dipole* payee=oldtow.front();
      if(payee->p_top==this->p_bot) payee->f_top=true; else payee->f_bot=true;
    }
    p_bot=new Antibranch(ati); assert(p_bot);
  }
  else p_bot=&ati;

  UpdateType();    //That concerns only this dipole.
  UpdateMass();

  list<Dipole*>& newtow=AccessBotBranch();
  newtow.remove(this);
  newtow.push_front(this);

}





void Dipole::RenewBranch(bool top, Glubranch& glu) {

  if(top) {    //Then method is concerned with p_top.

    if(p_top==&glu && !f_top) return;

    list<Dipole*>& oldtow=AccessTopBranch();
    oldtow.remove(this);

    if(f_top) {
      if(oldtow.empty()) { if(p_top) { delete p_top; p_top=NULL;}}
      else {
	Test(oldtow,this->p_top,string("Top"));
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

  }
  else {    //Then method is concerned with p_bot.

    if(p_bot==&glu && !f_bot) return;

    list<Dipole*>& oldtow=AccessBotBranch();
    oldtow.remove(this);

    if(f_bot) {
      if(oldtow.empty()) { if(p_bot) { delete p_bot; p_bot=NULL;}}
      else {
	Test(oldtow,this->p_bot,string("Bot"));
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





//eof
