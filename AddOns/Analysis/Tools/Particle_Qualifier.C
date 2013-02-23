#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Smart_Pointer.C"

using namespace ATOOLS;

namespace ATOOLS { template class SP(Particle_Qualifier_Base); }

Particle_Qualifier_Base::~Particle_Qualifier_Base()
{
}

void Particle_Qualifier_Base::Keep(Particle_List *const list)
{
  for (Particle_List::iterator pit=list->begin();pit!=list->end();) 
    if (!(*this)(*pit)) pit=list->erase(pit);
    else ++pit;
}

void Particle_Qualifier_Base::Erase(Particle_List *const list)
{
  for (Particle_List::iterator pit=list->begin();pit!=list->end();) 
    if ((*this)(*pit)) pit=list->erase(pit);
    else ++pit;
}

template void ATOOLS::copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_KF &);
template void ATOOLS::copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Final_State &);
template void ATOOLS::copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Charged &);

namespace ATOOLS {
  
  template<>
  Particle_Qualifier_Base *Getter_Function<Particle_Qualifier_Base,std::string>::
  GetObject(const std::string &name,const std::string &parameters)
  {
    msg_Tracking()<<"Looking for Qualifier '"<<name<<"' ... ";
    if (name[0]=='!') {
      std::string name1=name.substr(1);
      Particle_Qualifier_Base * qual = ATOOLS::Particle_Qualifier_Getter::GetObject(name1,name1);
      if (qual) return new Not_Particle_Qualifier(qual);
    }
    size_t pos=name.find("|");
    if (pos!=std::string::npos) {
      std::string name1=name.substr(0,pos);
      std::string name2=name.substr(pos+1);
      Particle_Qualifier_Base * qual1 = ATOOLS::Particle_Qualifier_Getter::GetObject(name1,name1);
      Particle_Qualifier_Base * qual2 = ATOOLS::Particle_Qualifier_Getter::GetObject(name2,name2);
      if (qual1 && qual2) return new Or_Particle_Qualifier(qual1,qual2);
    }
    pos=name.find("&");
    if (pos!=std::string::npos) {
      std::string name1=name.substr(0,pos);
      std::string name2=name.substr(pos+1);
      Particle_Qualifier_Base * qual1 = ATOOLS::Particle_Qualifier_Getter::GetObject(name1,name1);
      Particle_Qualifier_Base * qual2 = ATOOLS::Particle_Qualifier_Getter::GetObject(name2,name2);
      if (qual1 && qual2) return new And_Particle_Qualifier(qual1,qual2);
    }
    pos=name.find('(');
    if (pos!=std::string::npos) {
      size_t epos(name.rfind(')'));
      if (epos>pos) {
	return ATOOLS::Particle_Qualifier_Getter::
	  GetObject(name.substr(0,pos),name.substr(pos+1,epos-pos-1));
      }
    }
    String_Getter_Map::iterator git=s_getters->find(name);
    if (git!=s_getters->end()) {
      msg_Tracking()<<"found."<<std::endl;
      return (*git->second)(parameters);
    }
    msg_Tracking()<<"not found."<<std::endl;
    return NULL;
  }
}

#define COMPILE__Getter_Function
#define OBJECT_TYPE Particle_Qualifier_Base
#define PARAMETER_TYPE std::string
#include "ATOOLS/Org/Getter_Function.C"



template <class Class>
Particle_Qualifier_Base *GetQualifier(const std::string &parameter)
{									
  return new Class();
}									

#define DEFINE_GETTER_METHOD(CLASS)				\
  Particle_Qualifier_Base *					\
  ATOOLS::Getter<Particle_Qualifier_Base,std::string,CLASS>::	\
  operator()(const std::string &parameter) const		\
  { return GetQualifier<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void ATOOLS::Getter<Particle_Qualifier_Base,std::string,NAME>::	\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<PRINT; }

#define DEFINE_QUALIFIER_GETTER(CLASS,TAG,PRINT,DISP)		\
  DECLARE_ND_GETTER(CLASS,TAG,Particle_Qualifier_Base,std::string,DISP);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS,PRINT)

#include "ATOOLS/Org/Message.H"

void Particle_Qualifier_Base::ShowQualifiers(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<"Particle_Qualifier_Base::ShowQualifiers(): {\n\n";
  msg_Out()<<"   new qualifiers can be constructed\n";
  msg_Out()<<"   using the operators '!', '&' and '|'\n\n";
  Particle_Qualifier_Getter::PrintGetterInfo(msg->Out(),20);
  msg_Out()<<"\n}"<<std::endl;
}


#define DEFINE_QUALIFIER_CLASS(NAME)                              \
class NAME : public Particle_Qualifier_Base {                     \
public:                                                           \
  bool operator()(const Particle*) const;                         \
}


DEFINE_QUALIFIER_CLASS(Is_BHadron_Decay_Product);
DEFINE_QUALIFIER_GETTER(Is_BHadron_Decay_Product,
			"DecayedBHadron","decay product of bhadron",1)

DEFINE_QUALIFIER_CLASS(Is_BQuark_Decay_Product);
DEFINE_QUALIFIER_GETTER(Is_BQuark_Decay_Product,
			"DecayedBQuark","decayed b quark",1)

DEFINE_QUALIFIER_GETTER(Is_ME_Particle,
			"ME","ME particle",1)

DEFINE_QUALIFIER_GETTER(Is_Charged_Hadron,
			"ChargedHadron","charged hadron",1)
DEFINE_QUALIFIER_GETTER(Is_Neutral_Hadron,
			"NeutralHadron","neutral hadron",1)
DEFINE_QUALIFIER_GETTER(Is_Hadron,
			"Hadron","hadron",1)
DEFINE_QUALIFIER_GETTER(Is_Charged,
			"Charged","charged",1)
DEFINE_QUALIFIER_GETTER(Is_Neutral,
			"Neutral","neutral",1)
DEFINE_QUALIFIER_GETTER(Is_Charged_Pion,
			"ChargedPion","charged pion",0)
DEFINE_QUALIFIER_GETTER(Is_Charged_Kaon,
			"ChargedKaon","charged kaon",0)
DEFINE_QUALIFIER_GETTER(Is_Proton_Antiproton,
			"ProtonAntiproton","proton antiproton",0)
DEFINE_QUALIFIER_GETTER(Is_Parton,
			"Parton","parton",1)
DEFINE_QUALIFIER_GETTER(Is_There,
			"There","there",1)
DEFINE_QUALIFIER_GETTER(Is_Neutral_Pion,
			"NeutralPion","neutral pion",0)
DEFINE_QUALIFIER_GETTER(Is_Neutral_Kaon,
			"NeutralKaon","neutral kaon",0)
DEFINE_QUALIFIER_GETTER(Is_Charged_KStar,
			"ChargedKStar","charged kstar",0)
DEFINE_QUALIFIER_GETTER(Is_Neutral_KStar,
			"NeutralKStar","neutral kstar",0)
DEFINE_QUALIFIER_GETTER(Is_Eta,
			"Eta","eta",0)
DEFINE_QUALIFIER_GETTER(Is_Rho0,
			"Rho0","rho0",0)
DEFINE_QUALIFIER_GETTER(Is_Omega,
			"Omega","omega",0)
DEFINE_QUALIFIER_GETTER(Is_EtaPrime,
			"EtaPrime","eta prime",0)
DEFINE_QUALIFIER_GETTER(Is_Phi,
			"Phi","phi",0)
DEFINE_QUALIFIER_GETTER(Is_Lambda,
			"Lambda","lambda",0)
DEFINE_QUALIFIER_GETTER(Is_Charged_Sigma,
			"ChargedSigma","charged sigma",0)
DEFINE_QUALIFIER_GETTER(Is_Charged_Xi,
			"ChargedXi","charged xi",0)
DEFINE_QUALIFIER_GETTER(Is_Neutral_Xi,
			"NeutralXi","neutral xi",0)
DEFINE_QUALIFIER_GETTER(Is_Not_Lepton,
			"NotLepton","not lepton",1)

bool Or_Particle_Qualifier::operator() (const Particle * p) const {
  return ((*p_qual_a)(p) || (*p_qual_b)(p));
}
bool And_Particle_Qualifier::operator() (const Particle * p) const {
  return ((*p_qual_a)(p) && (*p_qual_b)(p));
}
bool Not_Particle_Qualifier::operator() (const Particle * p) const {
  return !(*p_qual_a)(p);
}

bool Is_ME_Particle::operator() (const Particle * p) const {
  if ( p && p->Info()=='H' ) return 1;
  return 0;  
}


bool Is_BHadron_Decay_Product::operator() (const Particle * p) const {
  if (!p) return 0;
  if (p->Flav().IsB_Hadron()) return 1;
  Blob * b = p->ProductionBlob();
  if (!b || b->NInP()!=1 || b->Type()==btp::Fragmentation) return 0;
  return operator()(b->InParticle(0));
}

bool Is_BQuark_Decay_Product::operator() (const Particle * p) const {
  if (!p) return 0;
  if (p->Flav().Kfcode()==kf_b) return 1;
  Blob * b = p->ProductionBlob();
  if (!b || b->Type()==btp::Beam || b->Type()==btp::Signal_Process) return 0;
  return operator()(b->InParticle(0));
}

DECLARE_GETTER(Is_KF,"KF",Particle_Qualifier_Base,std::string);
Particle_Qualifier_Base *ATOOLS::Getter<Particle_Qualifier_Base,std::string,Is_KF>::
operator()(const std::string &parameter) const  
{ return new Is_KF(parameter); }
void ATOOLS::Getter<Particle_Qualifier_Base,std::string,Is_KF>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"kf code, usage: KF(<kf code>)"; }

Is_KF::Is_KF(const std::string &kfcode):
  m_kfcode((kf_code)abs(ToType<int>(kfcode))) {}

bool Is_KF::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==m_kfcode ) return 1;
  return 0;
}

DECLARE_GETTER(Is_Flav,"Flav",Particle_Qualifier_Base,std::string);
Particle_Qualifier_Base *ATOOLS::Getter<Particle_Qualifier_Base,std::string,Is_Flav>::
operator()(const std::string &parameter) const  
{ return new Is_Flav(parameter); }
void ATOOLS::Getter<Particle_Qualifier_Base,std::string,Is_Flav>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"flavour, usage: Flav(<+- kf code>)"; }

Is_Flav::Is_Flav(const std::string &kfcode)
{
  int id(ToType<int>(kfcode));
  m_flav=Flavour((kf_code)abs(id)); 
  if (id<0) m_flav=m_flav.Bar();
}

bool Is_Flav::operator() (const Particle * p) const {
  if ( p && p->Flav()==m_flav ) return 1;
  return 0;
}

bool Is_Parton::operator() (const Particle * p)const {
  if ( p && ( p->Flav().IsGluon() || p->Flav().IsQuark() ) ) return 1;
  return 0;
}

bool Is_Charged::operator() (const Particle * p)const{
  if ( p && (p->Flav().IntCharge() !=0) ) return 1;
  return 0;
}

bool Is_Hadron::operator() (const Particle * p)const{
  if ( p && p->Flav().IsHadron() &&!p->Flav().IsDiQuark()) return 1;
  return 0;
}

bool Is_Charged_Lepton::operator() (const Particle * p)const{
  if ( p && p->Flav().IsLepton() && p->Flav().IntCharge()!=0) return 1;
  return 0;
}

bool Is_Charged_Hadron::operator() (const Particle * p)const{
  if ( p && p->Flav().IntCharge() !=0 &&
       p->Flav().IsHadron() &&!p->Flav().IsDiQuark()) return 1;
  return 0;
}

bool Is_Neutral_Hadron::operator() (const Particle * p)const{
  if ( p && p->Flav().IntCharge() ==0 &&
       p->Flav().IsHadron() &&!p->Flav().IsDiQuark()) return 1;
  return 0;
}

bool Is_Final_State::operator() (const Particle * p)const{
  if ( p && (p->Status() == 1) ) return 1;
  return 0;
}

bool Is_Neutral::operator() (const Particle * p)const{
  if ( p && p->Flav().IntCharge()==0) return 1;
  return 0;
}

bool Is_Charged_Pion::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_pi_plus) return 1;
  return 0;
}

bool Is_Neutral_Pion::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_pi) return 1;
  return 0;
}

bool Is_Charged_Kaon::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_K_plus) return 1;
  return 0;
}
bool Is_Neutral_Kaon::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_K) return 1;
  return 0;
}

bool Is_Charged_KStar::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_K_star_892_plus) return 1;
  return 0;
}

bool Is_Neutral_KStar::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_K_star_892) return 1;
  return 0;
}

bool Is_Rho0::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_rho_770) return 1;
  return 0;
}

bool Is_Eta::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_eta) return 1;
  return 0;
}

bool Is_EtaPrime::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_eta_prime_958) return 1;
  return 0;
}

bool Is_Phi::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_phi_1020) return 1;
  return 0;
}

bool Is_Omega::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_omega_782) return 1;
  return 0;
}

bool Is_Lambda::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_Lambda) return 1;
  return 0;
}

bool Is_Charged_Sigma::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_Sigma_minus) return 1;
  return 0;
}

bool Is_Charged_Xi::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_Xi_minus) return 1;
  return 0;
}

bool Is_Neutral_Xi::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==kf_Xi) return 1;
  return 0;
}

bool Is_Proton_Antiproton::operator() (const Particle * p) const {
  if ( p && p->Flav().Kfcode()==2212) return 1;
  return 0;
}

bool Is_Not_Lepton::operator() (const Particle * p) const {
  if ( p && !p->Flav().IsLepton() ) return 1;
  return 0;
}

bool Is_Not_Neutrino::operator() (const Particle * p) const {
  if ( p && !(p->Flav().IsLepton() && p->Flav().IntCharge()==0) ) return 1;
  return 0;
}

bool Is_There::operator() (const Particle * p) const {
  if ( p ) return 1;
  return 0;
}


class Is_Strong : public Particle_Qualifier_Base {
public:
  bool operator()(const Particle *) const;
};

bool Is_Strong::operator() (const Particle * p) const {
  if ( p && p->Flav().Strong() ) return 1;
  return 0;
}

DEFINE_QUALIFIER_GETTER(Is_Strong,
			"Strong","strong",1)
