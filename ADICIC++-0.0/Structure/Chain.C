//bof
//Version: 1 ADICIC++-0.0/2004/05/24

//Implementation of Chain.H.



#include <ios>
#include <iomanip>
#include <string>
#include <sstream>
#include "Chain.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



ostream& ADICIC::operator<<(ostream& ost, const ADICIC::Chain& cha) {
  unsigned i=0;
  ost<<endl<<"\e[1m"<<string(80,'=')<<"\e[0m"<<endl;
  for(list<Dipole*>::const_iterator it=cha.DipolePointerList().begin();
      it!=cha.DipolePointerList().end(); ++it) {
    ++i;
    ost<<setiosflags(ios::left)<<setw(4)<<"\e[33m\e[1m"<<i<<"\e[0m"
       <<resetiosflags(ios::left)<<(**it)<<endl;
  }
  ost<<"\e[1m"<<string(80,'=')<<"\e[0m"<<endl;
  ost<<"Chain type: ";
  switch(cha.ChainType()) {
  case 1: ost<<"Ring"; break;
  case 0: ost<<"Line"; break;
  default : ost<<"incorrect";
  }
  ost<<"    State: ";
  switch(cha.Status()) {
  case -1 : ost<<"blocked"; break;
  case 0  : ost<<"off"; break;
  default : ost<<"on";
  }
  ost<<"    Handler: "<<bool(cha.IsHandled())<<endl;
  ost<<setiosflags(ios::left)
     <<"Scale="<<setw(12)<<cha.LastScale()
     <<"Mass="<<setw(12)<<cha.Mass()
     <<"SqrMass="<<setw(12)<<cha.InvMass()
     <<"P="<<cha.Momentum();
  ost<<resetiosflags(ios::left)<<endl;
  ost<<"\e[1m"<<string(80,'=')<<"\e[0m";
  return ost;
}



//=============================================================================



const Chain::Initiator::Simple_Epem Chain::Initiator::simple_epem={};



//=============================================================================



int Chain::s_count=0;
const int& Chain::InStore=Chain::s_count;



//=============================================================================



Chain::Chain() {
  ++s_count;
  varset.Init();
}





Chain::Chain(const Chain& cha) {

  ++s_count;
  varset.Init();

  varset.f_active=cha.varset.f_active;
  varset.m_k2tlast=cha.varset.m_k2tlast;

  char key;

  Dipole::Branch* topb;
  Dipole::Antibranch* bota;
  Dipole::Glubranch* topg;
  Dipole::Glubranch* botg;

  for(list<Dipole*>::const_iterator dit=cha.varset.l_dip.begin();
      dit!=cha.varset.l_dip.end(); ++dit) {
    key=0; topb=NULL; bota=NULL; topg=botg=NULL;
    Dipole_Particle* to=(*dit)->GetTopBranchPointer().operator->();
    Dipole_Particle* bo=(*dit)->GetBotBranchPointer().operator->();
    if(to->OrgType()==Nil) topg=static_cast<Dipole::Glubranch*>(to);
    else                 { topb=static_cast<Dipole::Branch*>(to); key=2;}
    if(bo->OrgType()==Nil) botg=static_cast<Dipole::Glubranch*>(bo);
    else                 { bota=static_cast<Dipole::Antibranch*>(bo); key+=1;}
    if(topg) {
      topg=new Dipole::Glubranch(*topg);
      assert(topg);
      varset.l_glub.push_back(topg);
    } else {
      topb=new Dipole::Branch(*topb);
      assert(topb); assert(!varset.p_quab);
      varset.p_quab=topb;
    }
    if(botg) {
      botg=new Dipole::Glubranch(*botg);
      assert(botg);
      varset.l_glub.push_back(botg);
    } else {
      bota=new Dipole::Antibranch(*bota);
      assert(bota); assert(!varset.p_atib);
      varset.p_atib=bota;
    }
    Dipole* dip=NULL;
    switch(key) {
    case 0: dip=new Dipole(*topg,*botg); break;
    case 1: dip=new Dipole(*topg,*bota); break;
    case 2: dip=new Dipole(*topb,*botg); break;
    case 3: dip=new Dipole(*topb,*bota);
    }
    assert(dip);
    varset.l_dip.push_back(dip);
  }

  if(cha.CheckMomentum(varset.m_momentum)) {
    varset.m_mass=cha.varset.m_mass;
    varset.m_invmass=cha.varset.m_invmass;
  } else {
    varset.f_active=Blocked;
    varset.m_invmass=varset.m_momentum.Abs2();
    if(varset.m_invmass<0.0) {
      cerr<<"\nMethod: ADICIC::Chain::Chain(const ADICIC::Chain&): "
	  <<"Warning: Negative invariant mass!\n"<<endl;
      varset.m_mass=-1*sqrt(-1*varset.m_invmass);
    }
    else varset.m_mass=sqrt(varset.m_invmass);
  }

}





Chain::Chain(const Dipole::Branch& ban, const Dipole::Antibranch& ati,
	     const Initiator::Simple_Epem) {

  ++s_count;

  varset.p_hdl=NULL;
  varset.l_glub.clear();
  varset.l_dip.clear();

  varset.p_quab=new Dipole::Branch(ban);
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

}





Chain::Chain(const Dipole::Glubranch& glut, const Dipole::Glubranch& glub,
	     const Initiator::Simple_Epem) {

  ++s_count;

  if(&glut==&glub) {
    cerr<<"\nMethod: ADICIC::Chain::Chain(const ADICIC::Dipole::Glubranch&, "
	<<"const ADICIC::Dipole::Glubranch&, "
	<<"const ADICIC::Chain::Initiator::Simple_Epem): "
	<<"Warning: Gluon Chain initiation from a single gluon!\n"<<endl;
    //exit(1);
  }

  varset.p_hdl=NULL;
  varset.p_quab=NULL;
  varset.p_atib=NULL;
  varset.l_glub.clear();
  varset.l_dip.clear();

  Dipole::Glubranch* top=new Dipole::Glubranch(glut);
  assert(top);
  varset.l_glub.push_front(top);
  Dipole::Glubranch* bot=new Dipole::Glubranch(glub);
  assert(bot);
  varset.l_glub.push_back(bot);

  Dipole* root=new Dipole(*top,*bot);
  assert(root);
  varset.l_dip.push_front(root);

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

}





Chain::~Chain() {

  if(varset.p_hdl) {    //Check for Dipole_Handler.
    cerr<<"\nMethod: ADICIC::Chain::~Chain(): "
	<<"Warning: Detaching Chain_Handler from Dipole!\n"<<endl;
    //if(varset.p_hdl->IsDockedAt(*this)==false) {
    //  cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
    //  assert(varset.p_hdl->IsDockedAt(*this));
    //}
    Chain_Handler* phand=varset.p_hdl;
    varset.p_hdl=NULL;
    //phand->DetachDipole(this);
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

}



//=============================================================================



Chain& Chain::operator=(const Chain& cha) {

  if( this==&cha || !this->Clear() ) return *this;

  varset.f_active=cha.varset.f_active;
  varset.m_k2tlast=cha.varset.m_k2tlast;

  char key;

  Dipole::Branch*     topb;
  Dipole::Antibranch* bota;
  Dipole::Glubranch*  topg;
  Dipole::Glubranch*  botg;

  for(list<Dipole*>::const_iterator dit=cha.varset.l_dip.begin();
      dit!=cha.varset.l_dip.end(); ++dit) {
    key=0; topb=NULL; bota=NULL; topg=botg=NULL;
    Dipole_Particle* to=(*dit)->GetTopBranchPointer().operator->();
    Dipole_Particle* bo=(*dit)->GetBotBranchPointer().operator->();
    if(to->OrgType()==Nil) topg=static_cast<Dipole::Glubranch*>(to);
    else                 { topb=static_cast<Dipole::Branch*>(to); key=2;}
    if(bo->OrgType()==Nil) botg=static_cast<Dipole::Glubranch*>(bo);
    else                 { bota=static_cast<Dipole::Antibranch*>(bo); key+=1;}
    if(topg) {
      topg=new Dipole::Glubranch(*topg); assert(topg);
      varset.l_glub.push_back(topg);
    } else {
      topb=new Dipole::Branch(*topb); assert(topb); assert(!varset.p_quab);
      varset.p_quab=topb;
    }
    if(botg) {
      botg=new Dipole::Glubranch(*botg); assert(botg);
      varset.l_glub.push_back(botg);
    } else {
      bota=new Dipole::Antibranch(*bota); assert(bota); assert(!varset.p_atib);
      varset.p_atib=bota;
    }
    Dipole* dip=NULL;
    switch(key) {
    case 0: dip=new Dipole(*topg,*botg); break;
    case 1: dip=new Dipole(*topg,*bota); break;
    case 2: dip=new Dipole(*topb,*botg); break;
    case 3: dip=new Dipole(*topb,*bota);
    }
    assert(dip);
    varset.l_dip.push_back(dip);
  }

  if(cha.CheckMomentum(varset.m_momentum)) {
    varset.m_mass=cha.varset.m_mass;
    varset.m_invmass=cha.varset.m_invmass;
  } else {
    varset.f_active=Blocked;
    varset.m_invmass=varset.m_momentum.Abs2();
    if(varset.m_invmass<0.0) {
      cerr<<"\nMethod: ADICIC::Chain& ADICIC::Chain::operator=("
	  <<"const ADICIC::Chain&): "
	  <<"Warning: Negative invariant mass!\n"<<endl;
      varset.m_mass=-1*sqrt(-1*varset.m_invmass);
    }
    else varset.m_mass=sqrt(varset.m_invmass);
  }

  return *this;

}





//const bool Chain::operator==(const Chain& cha) const {
//  cerr<<"\nSorry :o( Method has not been implemented yet.\n";
//  assert(0); exit(1);
//}



//=============================================================================



void Chain::Print() const {
  cout<<(*this)<<endl;
  if(varset.p_quab) cout<<"Branch:"<<endl<<&(varset.p_quab->Parton)<<endl;
  if(varset.p_atib) cout<<"Antibranch:"<<endl<<&(varset.p_atib->Parton)<<endl;
  if(!varset.l_glub.empty()) cout<<"Glubranches:"<<endl;
  for(list<Dipole::Glubranch*>::const_iterator gut=varset.l_glub.begin();
      gut!=varset.l_glub.end(); ++gut) {
    cout<<&((*gut)->Parton)<<endl;
  }
  cout<<"\e[1m"<<string(80,'=')<<"\e[0m"<<endl;
  cout<<"Number of branches = "<<ParticleNumber()<<endl;
  cout<<"Number of  dipoles = "<<DipoleNumber()<<endl;
  cout<<"\e[1m"<<string(80,'=')<<"\e[0m"<<endl;
}





const bool Chain::CheckMomentum(ATOOLS::Vec4D& sum) const {
  sum=Vec4D();
  if(varset.p_quab) sum+=varset.p_quab->Momentum();
  if(varset.p_atib) sum+=varset.p_atib->Momentum();
  for(list<Dipole::Glubranch*>::const_iterator gut=varset.l_glub.begin();
      gut!=varset.l_glub.end(); ++gut)
    sum+=(*gut)->Momentum();
  if(sum==varset.m_momentum) return true;
  return false;
}





const bool Chain::Clear() {
  if(varset.p_hdl) {
    cerr<<"\nMethod: ADICIC::Chain::Clear(): Warning: "
	<<"Clearing a Chain manipulated by a Chain_Handler is not permitted!\n"
	<<endl;
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
  varset.Init();
  return true;
}





const bool Chain::Initialize(const Dipole::Branch& ban,
			     const Dipole::Antibranch& ati) {

  if(!Clear()) return false;

  varset.p_quab=new Dipole::Branch(ban);
  varset.p_atib=new Dipole::Antibranch(ati);
  assert(varset.p_quab);
  assert(varset.p_atib);

  Dipole* root=new Dipole(*varset.p_quab,*varset.p_atib);
  assert(root);
  varset.l_dip.push_front(root);

  varset.m_k2tlast=root->ProdScale();
  varset.m_mass=root->Mass();
  varset.m_invmass=root->InvMass();
  varset.m_momentum=root->TotP();

  varset.f_active=root->Status();
  if(varset.m_invmass<0.0 || varset.m_momentum[0]<0.0)
    varset.f_active=Blocked;

  return true;

}



//=============================================================================



const ATOOLS::Vec4D& Chain::UpdateMomentum(double k, const ATOOLS::Vec4D& p) {
  varset.m_momentum+=k*p;
  varset.m_invmass=varset.m_momentum.Abs2();
  if(varset.m_invmass<0.0) {
    cerr<<"\nMethod: const ATOOLS::Vec4D& ADICIC::Chain::UpdateMomentum("
	<<"double, const ATOOLS::Vec4D&): "
	<<"Warning: Negative invariant mass!\n"<<endl;
    varset.m_mass=-1*sqrt(-1*varset.m_invmass);
  }
  else varset.m_mass=sqrt(varset.m_invmass);
  return varset.m_momentum;
}



//=============================================================================





//eof
