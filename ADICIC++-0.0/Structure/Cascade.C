//bof
//Version: 2 ADICIC++-0.0/2005/01/31

//Implementation of Cascade.H.



#ifdef __GNUC__
#if __GNUC__ >2
#include <ios>
#endif
#endif
#include <iomanip>
#include <string>
#include <sstream>
#include "Message.H"
#include "Cascade.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



ostream& ADICIC::operator<<(ostream& ost, const ADICIC::Cascade& cas) {
  ost<<"\n"<<om::bold<<string(140,'O')<<om::reset<<"\n";
  ost<<"Cascade type: ";
  switch(cas.CascadeType()) {
  case 10: ost<<"lines"; break;
  case 0 : ost<<"line"; break;
  case 1 : ost<<"ring"; break;
  default: ost<<"incorrect";
  }
  ost<<"    Cascade state: ";
  switch(cas.Status()) {
  case -1: ost<<"blocked"; break;
  case 0 : ost<<"off"; break;
  default: ost<<"on";
  }
  ost<<"    Cascade handler: "<<cas.IsHandled();
  ost<<"    New chain in preparation: "<<cas.ChainPreparationPointer()<<"\n";
  ost<<setiosflags(ios::left)
     <<"Mass="<<setw(15)<<cas.Mass()
     <<"SqrMass="<<setw(15)<<cas.InvMass()
     <<"Momentum="<<cas.Momentum();
  ost<<resetiosflags(ios::left)<<"\n";
  if(Cascade::s_print) {
    ost<<om::bold<<string(140,'~')<<om::reset<<"\n";
    ost<<"Number of root chains       = "<<cas.RootChainNumber()<<"\n";
    ost<<"Number of chains in cascade = "<<cas.ChainNumber()<<"\n";
    ost<<"Total number of branches    = "<<cas.ParticleNumber()<<"\n";
    ost<<"Total number of  dipoles    = "<<cas.DipoleNumber()<<"\n";
  }
  ost<<om::bold<<string(140,'~')<<om::reset<<"\n";
  list<Cascade::Mirror>::const_iterator mit=cas.MirrorList().begin();
  for(list<Chain*>::const_iterator cit=cas.ChainPointerList().begin();
      cit!=cas.ChainPointerList().end(); ++cit) {
    if((*mit).first==0) ost<<om::bold<<om::red;
    else ost<<om::bold<<om::blue;
    ost<<setiosflags(ios::left)<<"Mother chain = "
       <<setw(10)<<(*mit).first<<"Number of splittings of this chain = "
       <<setw(10)<<(*mit).second
       <<om::reset<<resetiosflags(ios::left);
    if(Cascade::s_print) {
      bool temp=Chain::MoreCout;
      Chain::MoreCout=true;
      ost<<(**cit)<<endl;
      Chain::MoreCout=temp;
    }
    else ost<<(**cit)<<endl;
    ++mit;
  }
  ost<<om::bold<<string(140,'O')<<om::reset;
  return ost;
}



//=============================================================================



//const Cascade::Initiator::Simple_EpEm Cascade::Initiator::simple_epem={};



//=============================================================================



bool Cascade::s_print=false;

int Cascade::s_count=0;
const int& Cascade::InStore=Cascade::s_count;



//=============================================================================



Cascade::Cascade() {
  ++s_count;
  caset.Init();    //This is a clear cascade.
}





Cascade::Cascade(const Cascade& cas) {
  ++s_count;
  caset.Init();
  //Later: Additional Chain and Cascade initialization methods will cause
  //further modifications here!
  Type type=cas.CascadeType();
  if(type==incorrect) return;//////////////////////////////////////////////////
  this->Copy(cas,type);
}





Cascade::~Cascade() {

  if(caset.p_hdl) {    //Check for Cascade_Handler.
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Detaching Cascade_Handler from Cascade!\n"<<endl;
    if(caset.p_hdl->IsDockedAt(*this)==false) {
      cerr<<"\nBug: Wrong Cascade-Cascade_Handler connection emerged!\n";
      assert(caset.p_hdl->IsDockedAt(*this));
    }
    Cascade_Handler* phand=caset.p_hdl;
    caset.p_hdl=NULL;
    phand->DetachCascade(this);
  }

  if(caset.p_add) delete caset.p_add;
  for(list<Chain*>::iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit)
    if(*cit) delete (*cit);

  caset.l_ifo.clear();
  caset.l_cha.clear();

  --s_count;

}



//=============================================================================



Cascade& Cascade::operator=(const Cascade& cas) {
  if(this==&cas) return *this;
  //Later: Additional Chain and Cascade initialization methods will cause
  //further modifications here!
  Type type=cas.CascadeType();
  if(type==incorrect) return *this;////////////////////////////////////////////
  if(this->Clear()==false) return *this;
  this->Copy(cas,type);
  return *this;
}





//const bool Cascade::operator==(const Cascade& cha) const {
//  cerr<<"\nSorry :o( Method has not been implemented yet.\n";
//  assert(0); exit(1);
//}



//=============================================================================



const Cascade::Type Cascade::CascadeType() const {
  size_t num=caset.l_cha.size();
  if(num==0) return Cascade::incorrect;
  if(num==1) {
    Chain::Type typ=caset.l_cha.front()->ChainType();
    if(typ==Chain::line) return Cascade::line;
    if(typ==Chain::ring) return Cascade::ring;
    return Cascade::incorrect;
  }
  for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
      cat!=caset.l_cha.end(); ++cat)
    if((*cat)->ChainType()!=Chain::line) return Cascade::incorrect;
  return Cascade::lines;
}





const bool Cascade::CheckMomentumConservation(ATOOLS::Vec4D& sum) const {
  sum=Vec4D();
  Vec4D temp;
  for(list<Chain*>::const_iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit) {
    (*cit)->CheckMomentumConservation(temp);
    sum+=temp;
  }
  if(sum==caset.m_momentum) return true;
  return false;
}





const bool Cascade::ExtractPartons(list<Particle_List>& lists) const {

  Type type=this->CascadeType();
  if(type==incorrect) return false;

  if(lists.empty());
  else {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
	<<"List of particle pointer lists has already an entry!\n"<<endl;
  }

  size_t root=0;
  list<Particle_List>::iterator lit;
  list<Chain*>::const_iterator cit;
  for(cit=caset.l_cha.begin(); cit!=caset.l_cha.end(); ++cit) {
    const Chain& cha=**cit;
    if(cha.Memo!=root) {
      root=cha.Memo;
      lit=lists.insert(lists.end(),Particle_List());
    }
    assert(cha.ExtractPartons(*lit));
  }

#ifdef CASCADE_OUTPUT
  for(list<Particle_List>::const_iterator loc=lists.begin();
      loc!=lists.end(); ++loc)
    cout<<(*loc)<<endl;
#endif

  return true;

}





const bool Cascade::Clear() {
  if(caset.p_hdl) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Clearing a Cascade "
	<<"manipulated by a Cascade_Handler is not permitted!\n"<<endl;
    return false;
  }
  if(caset.p_add) delete caset.p_add;
  for(list<Chain*>::iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit)
    if(*cit) delete (*cit);
  caset.Init();    //This is a clear cascade right now.
  return true;
}





const bool Cascade::AddChain(const Dipole::Branch& ban,
			     const Dipole::Antibranch& ati,
			     bool ontop) {
  if(caset.p_add) return false;

  caset.p_add=new Chain(); assert(caset.p_add);
  if(caset.p_add->Initialize(ban,ati)==false) {
    delete caset.p_add;
    caset.p_add=NULL;
    return false;
  }

  if(ontop) {
    caset.l_ifo.push_front(Mirror(0,0));
    caset.l_cha.push_front(caset.p_add);
  } else {
    caset.l_ifo.push_back(Mirror(0,0));
    caset.l_cha.push_back(caset.p_add);
  }

  UpdateMomentum(1.0,caset.p_add->Momentum());

  if(caset.m_nroot) {
    if(caset.f_active==On) caset.f_active=caset.p_add->Status();
  } else {
    caset.f_active=caset.p_add->Status();
  }
  if(caset.m_invmass<0.0 || caset.m_momentum[0]<0.0) caset.f_active=Blocked;

  ++caset.m_nroot;
  caset.p_add=NULL;

  return true;
}




/*
const bool Cascade::AddChain(const Dipole::Glubranch& glut,
			     const Dipole::Glubranch& glub) {

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
*/


//=============================================================================



void Cascade::Copy(const Cascade& cac, const Type type) {

  //Later: Additional Chain and Cascade initialization methods will cause
  //further modifications here!

  //Only two types are now possible: line and ring.
  //A line has one dipole at least.
  //A ring has two dipoles at least and the root gluon always resides at the
  //first place of the Glubranch list.

  //The variable type is not needed yet.

  caset.f_active=cac.caset.f_active;
  caset.m_nroot=cac.caset.m_nroot;

  caset.l_ifo=cac.caset.l_ifo;
  caset.l_cha=cac.caset.l_cha;
  //caset.l_ifo.assign(cac.caset.l_ifo.begin(),cac.caset.l_ifo.end());
  //caset.l_cha.assign(cac.caset.l_cha.begin(),cac.caset.l_cha.end());

  ////////////////////////////////////testing was//////////////////////////////

  //Real copies and mother chain name correction.
  map<size_t,size_t> link; link[0]=0;
  list<Chain*>::iterator cit=caset.l_cha.begin();
  for(; cit!=caset.l_cha.end(); ++cit) {
    size_t oldname=(*cit)->Name;
    *cit=new Chain(**cit); assert(*cit);
    link[oldname]=(*cit)->Name;
  }
  list<Mirror>::iterator mit=caset.l_ifo.begin();
  for(; mit!=caset.l_ifo.end(); ++mit) (*mit).first=link[(*mit).first];

  /////////////////////////////////////////////////////////////////////////////

  if(cac.CheckMomentumConservation(caset.m_momentum)) {
    caset.m_mass=cac.caset.m_mass;
    caset.m_invmass=cac.caset.m_invmass;
  } else {
    caset.f_active=Blocked;
    if(caset.m_momentum[0]<0.0)
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Warning: Total energy is negative!\n"<<endl;
    caset.m_invmass=caset.m_momentum.Abs2();
    if(caset.m_invmass<0.0) {
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Warning: Negative invariant mass!\n"<<endl;
      caset.m_mass=-1*sqrt(-1*caset.m_invmass);
    }
    else caset.m_mass=sqrt(caset.m_invmass);
  }

}





const ATOOLS::Vec4D&
Cascade::UpdateMomentum(double k, const ATOOLS::Vec4D& p) {
  caset.m_momentum+=k*p;
  caset.m_invmass=caset.m_momentum.Abs2();
  if(caset.m_invmass<-1.0e-10) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
	<<"Negative invariant mass ("<<caset.m_invmass<<") !\n"<<endl;
    //caset.f_active=Blocked;
    caset.m_mass=-1*sqrt(-1*caset.m_invmass);
  }
  else caset.m_mass=sqrt(caset.m_invmass);
  return caset.m_momentum;
}



//=============================================================================





//eof
