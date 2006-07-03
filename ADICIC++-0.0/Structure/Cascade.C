//bof
//Version: 4 ADICIC++-0.0/2006/05/23

//Implementation of Cascade.H.



#ifdef __GNUC__
#if __GNUC__ >2
#include <ios>
#endif
#endif
#include <iostream>
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
  int pcs=ost.precision(6);
  ost<<"\n"<<om::bold<<string(130,'O')<<om::reset<<"\n";
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
    ost<<om::bold<<string(130,'~')<<om::reset<<"\n";
    ost<<setiosflags(ios::left);
    ost<<"Number of      : root chains = "<<setw(10)<<cas.RootChainNumber();
    ost<<"chains   = "<<setw(10)<<cas.ChainNumber()<<"\n";
    ost<<"Total number of: dipoles     = "<<setw(10)<<cas.DipoleNumber();
    ost<<"branches = "<<setw(10)<<cas.ParticleNumber();
    ost<<"related particles = "<<setw(10)<<cas.RelatedParticleNumber()<<"\n";
    ost<<resetiosflags(ios::left);
  }
  ost<<om::bold<<string(130,'~')<<om::reset<<"\n";
  list<Cascade::Mirror>::const_iterator mit=cas.MirrorList().begin();
  for(list<Chain*>::const_iterator cit=cas.ChainPointerList().begin();
      cit!=cas.ChainPointerList().end(); ++cit) {
    if((*mit).first==(*cit)->Name) ost<<om::bold<<om::red;
    else ost<<om::bold<<om::blue;
    ost<<setiosflags(ios::left)<<"Root chain = "
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
  ost<<om::bold<<string(130,'O')<<om::reset;
  ost.precision(pcs);
  return ost;
}



//=============================================================================



//const Cascade::Initiator::EpEm Cascade::Initiator::epem={};



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
  if(type==incorrect) return;
  this->Copy(cas);
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
  this->Destruct();
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
  //The message will come from Clear() itself.
  this->Copy(cas);
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





const unsigned Cascade::INumber() const {
  unsigned in=0;
  for(list<Chain*>::const_iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit)
    in+=(*cit)->INumber();
  return in;
}





const bool Cascade::CheckMomentumConservation(Vec4D& sum) const {
  sum=Vec4D();
  Vec4D test;
  for(list<Chain*>::const_iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit) {
    (*cit)->CheckMomentumConservation(test);
    sum+=test;
  }
  if(sum==caset.m_momentum) return true;
  test=sum-caset.m_momentum;
  for(char i=0; i<4; ++i) if(dabs(test[i])>1.0e-10) return false;
  return true;
}





const boolint Cascade::ExtractPartons(list<Particle_List>& lists) const {

  Type type=this->CascadeType();
  if(type==incorrect) return boolint();    //(false,0)!

  if(lists.empty());
  else {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
	<<"List of particle pointer lists has already an entry!\n"<<endl;
  }

  //Assume evolutions are separated in subsequent root-chain blocks.
  //One particle list per root chain including its correlated particles
  //plus its daughters.
  size_t root=0;
  int num=32767;
  list<Mirror>::const_iterator mit=caset.l_ifo.begin();
  list<Chain*>::const_iterator cit=caset.l_cha.begin();
  list<Particle_List>::iterator lit;
  for(; cit!=caset.l_cha.end(); ++cit) {
    const Chain& cha=**cit;
    size_t aroot=(*mit).first; assert(aroot!=0);
    if(aroot!=root) {    //In the first time this is always fulfilled.
      root=aroot;
      lit=lists.insert(lists.end(),Particle_List());
    }
    boolint res=cha.ExtractPartons(*lit);
    assert(res.flag);
    if(res.numr<num) num=res.numr;
    ++mit;
  }

#ifdef CASCADE_OUTPUT
  for(list<Particle_List>::const_iterator loc=lists.begin();
      loc!=lists.end(); ++loc)
    cout<<(*loc)<<endl;
#endif

  return boolint(true,num);

}





const bool Cascade::Clear() {
  if(caset.p_hdl) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: Clearing a Cascade "
	<<"manipulated by a Cascade_Handler is not permitted!\n"<<endl;
    return false;
  }
  this->Destruct();
  caset.Init();    //This is a clear cascade right now.
  return true;
}





const bool Cascade::AddChain(const Dipole::Branch& ban,
			     const Dipole::Antibranch& ati,
			     const Particle_List* pleps,
			     double scale, bool ontop) {
  if(caset.p_add) return false;

  caset.p_add=new Chain(); assert(caset.p_add);
  if(pleps) {
    for(Particle_List::const_iterator pit=pleps->begin();
	pit!=pleps->end(); ++pit)
      assert(caset.p_add->AddCorrParticleToInit(**pit));
  }
  if(caset.p_add->Initialize(ban,ati,scale)==false) {
    delete caset.p_add;
    caset.p_add=NULL;
    return false;
  }

  if(ontop) {
    caset.l_ifo.push_front(Mirror(caset.p_add->Name,0));
    caset.l_cha.push_front(caset.p_add);
  } else {
    caset.l_ifo.push_back(Mirror(caset.p_add->Name,0));
    caset.l_cha.push_back(caset.p_add);
  }

  UpdateMomentum(1.0,caset.p_add->Momentum());

  if(caset.m_nroot) {
    if(caset.f_active==On) caset.f_active=caset.p_add->Status();
  } else {
    caset.f_active=caset.p_add->Status();
  }

  //if(caset.m_invmass<0.0||caset.m_momentum[0]<0.0) caset.f_active=Blocked;???

  ++caset.m_nroot;
  caset.p_add=NULL;

  return true;
}




/*
const bool Cascade::AddChain(const Dipole::Glubranch& glut,
                             const Dipole::Glubranch& glub,
                             const Particle_List* pleps,
                             double scale, bool ontop) {

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



void Cascade::Destruct() {
  if(caset.p_add) delete caset.p_add;
  for(std::list<Chain*>::iterator cit=caset.l_cha.begin();
      cit!=caset.l_cha.end(); ++cit)
    if(*cit) delete (*cit);
  caset.l_ifo.clear(); caset.l_cha.clear();
}





void Cascade::Copy(const Cascade& cac) {

  //Later: Additional Cascade initialization methods might cause further
  //modifications here!

  caset.f_active=cac.caset.f_active;
  caset.m_nroot=cac.caset.m_nroot;

  caset.l_ifo=cac.caset.l_ifo;
  caset.l_cha=cac.caset.l_cha;    //Only pointers are copied.
  //caset.l_ifo.assign(cac.caset.l_ifo.begin(),cac.caset.l_ifo.end());
  //caset.l_cha.assign(cac.caset.l_cha.begin(),cac.caset.l_cha.end());

  //Real copy creation and root/source chain name corrections.
  map<size_t,size_t> link; link[0]=0;
  list<Chain*>::iterator cit=caset.l_cha.begin();
  for(; cit!=caset.l_cha.end(); ++cit) {
    size_t oldname=(*cit)->Name;
    *cit=new Chain(**cit); assert(*cit);
    link[oldname]=(*cit)->Name;
  }
  for(cit=caset.l_cha.begin(); cit!=caset.l_cha.end(); ++cit)
    (*cit)->SetSource()=link[(*cit)->Source];
  list<Mirror>::iterator mit=caset.l_ifo.begin();
  for(; mit!=caset.l_ifo.end(); ++mit) (*mit).first=link[(*mit).first];

  //Test 4-momenta.
  if(cac.CheckMomentumConservation(caset.m_momentum)) {
    caset.m_mass=cac.caset.m_mass;
    caset.m_invmass=cac.caset.m_invmass;
  } else {
    //Probably, this part is not really needed.
    caset.f_active=Blocked;
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Momentum conservation check failed!\n"<<endl;
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





const Vec4D& Cascade::UpdateMomentum(double k, const Vec4D& p) {
  caset.m_momentum+=k*p;
  caset.m_invmass=caset.m_momentum.Abs2();
  if(caset.m_invmass<0.0) {
    //cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": Warning: "
    //    <<"Negative invariant mass ("<<caset.m_invmass<<") !\n"<<endl;
    //caset.f_active=Blocked;//????????????????????????????????????????????????
    caset.m_mass=-1*sqrt(-1*caset.m_invmass);
    //The minus sign functions only as a flag.
  }
  else caset.m_mass=sqrt(caset.m_invmass);
  return caset.m_momentum;
}



//=============================================================================





//eof
