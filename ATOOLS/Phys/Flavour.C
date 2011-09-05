#include "ATOOLS/Phys/Flavour.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

namespace ATOOLS 
{
  KF_Table s_kftable;
}

using namespace ATOOLS;

Particle_Info::Particle_Info(const Particle_Info &info):
  m_kfc(info.m_kfc), m_mass(info.m_mass), m_hmass(info.m_hmass),
  m_yuk(info.m_yuk), m_width(info.m_width),
  m_dg(info.m_dg), m_dm(info.m_dm), m_qoverp2(info.m_qoverp2), 
  m_icharge(info.m_icharge), m_isoweak(info.m_isoweak), 
  m_strong(info.m_strong), m_resummed(info.m_resummed),
  m_spin(info.m_spin), m_stable(info.m_stable), 
  m_masssign(info.m_masssign), m_dummy(info.m_dummy), m_majorana(info.m_majorana), 
  m_formfactor(0), m_on(info.m_on), m_massive(info.m_massive), m_hadron(info.m_hadron),
  m_idname(info.m_idname), m_texname(info.m_texname)
{
  m_content.resize(info.m_content.size());
  for (size_t i(0);i<info.m_content.size();++i) 
    m_content[i] = new Flavour(*info.m_content[i]); 
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass,const double &width,
 const int icharge,const int isoweak,const int strong,
 const int spin,const int majorana,const bool on,
 const int stable,bool massive,const std::string &idname,
 const std::string &texname,const bool dummy):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_yuk(mass), m_width(width),
  m_dg(0.0), m_dm(0.0), m_qoverp2(1.0), m_icharge(icharge),
  m_isoweak(isoweak), m_strong(strong), m_resummed(0), m_spin(spin), 
  m_stable(stable), m_masssign(1), m_dummy(dummy), m_majorana(majorana), 
  m_formfactor(0), m_on(on), m_massive(massive), m_hadron(0), m_idname(idname),
  m_texname(texname)
{
  m_content.push_back(new Flavour(*this));
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass,const double &width,
 const int icharge,const int isoweak,const int spin,const bool on,
 const int stable,const std::string &idname,const std::string &texname):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_yuk(0.0), m_width(width),
  m_icharge(icharge), m_isoweak(isoweak), m_strong(0), m_resummed(0), m_spin(spin), 
  m_stable(stable), m_masssign(1), m_dummy(0), m_majorana(0), m_formfactor(0), m_on(on),
  m_massive(1), m_hadron(1), m_idname(idname), m_texname(texname) 
{
  m_content.push_back(new Flavour(*this));
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass, const int icharge, const int spin,
 const int formfactor, const std::string &idname, const std::string &texname):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_yuk(0.0), m_width(0),
  m_icharge(icharge), m_isoweak(0), m_strong(0), m_resummed(0), m_spin(0),
  m_stable(1), m_masssign(1), m_dummy(0), m_majorana(0),  m_formfactor(formfactor),
  m_on(1), m_massive(1), m_hadron(1), m_idname(idname), m_texname(texname)
{
  m_content.push_back(new Flavour(*this));
}

Particle_Info::~Particle_Info()
{ 
  Clear();
}

void Particle_Info::Clear()
{ 
  for (size_t i(0);i<m_content.size();++i) 
    delete m_content[i]; 
  m_content.clear();
}

void Particle_Info::SetResummed()
{
  for (size_t i(0);i<m_content.size();++i) {
    s_kftable[m_content[i]->Kfcode()]->m_resummed=true;
    s_kftable[kf_resummed]->Add(*m_content[i]);
  }
}

void Particle_Info::Add(const Flavour &fl)
{ 
  for (size_t i(0);i<fl.Size();++i) 
    m_content.push_back(new Flavour(fl[i])); 
}

Flavour Particle_Info::operator[](const size_t &i) const
{ 
  return *m_content[i]; 
}

bool Particle_Info::Includes(const Flavour &fl) const
{
  for (size_t j(0);j<fl.Size();++j) {
    bool found(false);
    for (size_t i(0);i<m_content.size();++i)
      if (m_content[i]->Kfcode()==fl[j].Kfcode()) { 
	found=true; 
	break; 
      }
    if (!found) return false;
  }
  return true;
}

KF_Table::~KF_Table()
{
  for (const_iterator kfit(begin());kfit!=end();++kfit)
    delete kfit->second;
}

kf_code KF_Table::KFFromIDName(const std::string &idname) const
{
  for(const_iterator kfit(begin());kfit!=end();++kfit) 
    if (kfit->second->m_idname==idname) return kfit->first;
  return kf_none;
}

kf_code KF_Table::KFFromTexName(const std::string &texname) const
{
  for(const_iterator kfit(begin());
      kfit!=end();++kfit) 
    if (kfit->second->m_texname==texname) return kfit->first;
  return kf_none;
}

int Flavour::Ctq() const
{
  long int code(*this);
  if (IsGluon()) return 0;
  return code+6;
}

void Flavour::FromCtq(const int code)
{
  m_anti=code<6;
  if (code==6) p_info=s_kftable[kf_gluon];
  else p_info=s_kftable[(kf_code)abs(code-6)];
}

int Flavour::HepEvt() 
{
  switch (Kfcode()) {
  case kf_a_0_1450:      return 10111;
  case kf_a_0_1450_plus: return m_anti?-10211:10211;
  case kf_f_0_1370:      return 10221;
  case kf_f_0_1710:      return 10331;
  case kf_a_0_980:       return 9000111;
  case kf_a_0_980_plus:  return m_anti?-9000211:9000211;
  case kf_f_0_980:       return 9010221;
  case 13122:            return 23122;
  case 23122:            return 13122;
  }
  return (long int)*this;
}

void Flavour::FromHepEvt(long int code) 
{
  m_anti=code<0;
  code=(kf_code)abs(code);
  p_info=s_kftable[PdgToSherpa(code)];
}

std::string Flavour::TexName() const 
{
  std::string name;
  if (m_anti && (!SelfAnti())) name="\\overline{";
  switch (Kfcode()) {
  case kf_pi : {name+=std::string("\\pi^{0}");break;}
  case kf_K : {name+=std::string("K^{0}");break;}
  case kf_K_L : {name+=std::string("K_{L}");break;}
  case kf_K_S : {name+=std::string("K_{S}");break;}
  case kf_D : {name+=std::string("D^{0}"); break;}
  case kf_B : {name+=std::string("B^{0}"); break;}
  default :
    if (IsHadron()) {
      name+=IDName();
      name=StringReplace(name, "+", "{^{+}}");
      name=StringReplace(name, "-", "{^{-}}");
      name=StringReplace(name, "*", "{^{*}}");
      name=StringReplace(name, "'", "{^{\\prime}}");
      name=StringReplace(name, "(", "{_{(");
      name=StringReplace(name, ")", ")}}");
      name=StringReplace(name, "eta", "\\eta ");
      name=StringReplace(name, "rho", "\\rho ");
      name=StringReplace(name, "pi", "\\pi ");
      name=StringReplace(name, "omega", "\\omega ");
      name=StringReplace(name, "phi", "\\phi ");
      name=StringReplace(name, "psi", "\\psi ");
      name=StringReplace(name, "Delta", "\\Delta ");
      name=StringReplace(name, "Sigma", "\\Sigma ");
      name=StringReplace(name, "Lambda", "\\Lambda ");
      name=StringReplace(name, "Xi", "\\Xi ");
      name=StringReplace(name, "Omega", "\\Omega ");
    }
    else {
      name+=p_info->m_texname;
    }
    break;
  }
  if (m_anti && (!SelfAnti())) {
    name+="}";
    switch (Kfcode()) {
    case kf_pi : {name=std::string("\\pi^{0}");break;}
    case kf_pi_plus : {name=std::string("\\pi^{-}");break;}
    case kf_K : {name=std::string("\\bar K^{0}");break;}
    case kf_K_plus : {name=std::string("K^{-}");break;}
    case kf_K_L : {name=std::string("K_{L}");break;}
    case kf_K_S : {name=std::string("K_{S}");break;}
    case kf_D_plus : {name=std::string("D^{-}"); break;}
    case kf_B_plus : {name=std::string("B^{-}"); break;}
    case kf_rho_770_plus : {name=std::string("\\rho_{(770)}^{-}"); break;}
    case kf_rho_1450_plus : {name=std::string("\\rho_{(1450)}^{-}"); break;}
    case kf_rho_1700_plus : {name=std::string("\\rho_{(1700)}^{-}"); break;}
    }
  }
  return name;
}

std::string Flavour::RootName() const
{
  std::string name=StringReplace(TexName(), "\\", "#");
  name=StringReplace(name,"#overline","#bar");
  return name;
}

std::string Flavour::ShellName() const 
{
  std::string name(IDName());
  size_t pos(0);
  while ((pos=name.find("("))!=std::string::npos) name.replace(pos,1,"");
  while ((pos=name.find(")"))!=std::string::npos) name.replace(pos,1,"");
  while ((pos=name.find("/"))!=std::string::npos) name.replace(pos,1,"");
  while ((pos=name.find("'"))!=std::string::npos) name.replace(pos,1,"prime");
  while ((pos=name.find("*"))!=std::string::npos) name.replace(pos,1,"star");
  return name;
}

std::string Flavour::IDName() const 
{
  std::string name(p_info->m_idname);
  if (Kfcode()==kf_e || Kfcode()==kf_mu || Kfcode()==kf_tau) {
    name.erase(name.length()-1,1);
    if (IsAnti()) name+="+";
    else name+="-";      
  }
  else {
    if (Kfcode()==kf_Hplus || Kfcode()==kf_Wplus ||
        (IsHadron() && !IsBaryon() && name[name.length()-1]=='+')) {
      name.erase(name.length()-1,1);
      if (IsAnti()) name+="-";
      else name+="+";      
    }
    else if (IsAnti()) name+="b"; 
  }
  return name;
}

bool Flavour::IsDiQuark() const 
{
  if(abs(Kfcode())>=1103&&abs(Kfcode())<=5505) {
    double help=abs(Kfcode())/100.0-int(abs(Kfcode())/100.0); 
    if(help<0.031) return true;
  }
  return false;
}

bool Flavour::IsBaryon() const 
{
  if (abs(Kfcode())%10000<1000) return false;
  return !IsDiQuark();
}

bool Flavour::IsB_Hadron() const 
{
  if (abs(Kfcode())<100)                            return 0;
  if (Kfcode()-100*int(Kfcode()/100)<10)                 return 0;
  if (abs((Kfcode()-100*int(Kfcode()/100))/10)==5)       return 1;
  if (abs((Kfcode()-1000*int(Kfcode()/1000))/100)==5)    return 1;
  if (abs((Kfcode()-10000*int(Kfcode()/10000))/1000)==5) return 1;
  return 0;
}

bool Flavour::IsC_Hadron() const 
{
  if (abs(Kfcode())<100)                            return 0;
  if (Kfcode()-100*int(Kfcode()/100)<10)                 return 0;
  if (abs((Kfcode()-100*int(Kfcode()/100))/10)==4)       return 1;
  if (abs((Kfcode()-1000*int(Kfcode()/1000))/100)==4)    return 1;
  if (abs((Kfcode()-10000*int(Kfcode()/10000))/1000)==4) return 1;
  return 0;
}

double Flavour::GenerateLifeTime() const
{
  double proper_time = rpa->hBar() / Width();
  return -proper_time*log(1.-ran->Get());
}

bool Flavour::IsStable() const
{
  if (p_info->m_stable==0) return false;
  if (p_info->m_stable==1) return true;
  if (p_info->m_stable==2 && !IsAnti()) return true;
  if (p_info->m_stable==3 && IsAnti()) return true;
  return false;
}

kf_code Flavour::PdgToSherpa(const unsigned long& pdg)
{
  switch (pdg) {
  case 10111:   return kf_a_0_1450; break;
  case 10211:   return kf_a_0_1450_plus; break;
  case 10221:   return kf_f_0_1370; break;
  case 10331:   return kf_f_0_1710; break;
  case 13122:   return 23122; break;
  case 23122:   return 13122; break;
  case 9000111: return kf_a_0_980; break;
  case 9000211: return kf_a_0_980_plus; break;
  case 9010221: return kf_f_0_980; break; 
  case 91:      return kf_cluster; break;
  case 92:      return kf_string; break;
  default:      return pdg; break;
  }
}

std::ostream &ATOOLS::operator<<(std::ostream &os,const Flavour &fl)
{
  return os<<fl.IDName();
}

void ATOOLS::OutputHadrons(std::ostream &str) {
  
  str<<" List of Hadron data \n";
  str<<"          IDName";
  str<<std::setw(10)<<"kfc";
  str<<std::setw(14)<<"MASS[<kfc>]";
  str<<std::setw(16)<<"WIDTH[<kfc>]";
  str<<std::setw(16)<<"STABLE[<kfc>]";
  str<<std::setw(17)<<"ACTIVE[<kfc>]\n";
  
  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();
  
  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if ((flav.IsHadron() || flav.IsDiQuark())
	&& flav.Size()==1 && flav.Kfcode()!=0) {
      str<<std::setw(16)<<flav.IDName();
      str<<std::setw(10)<<flav.Kfcode();
      str<<std::setw(14)<<flav.HadMass();
      str<<std::setw(16)<<flav.Width();
      str<<std::setw(16)<<flav.Stable();
      str<<std::setw(16)<<flav.IsOn();
      str<<"\n";
    }
  }  
  str<<"\n";
}
     
void ATOOLS::OutputParticles(std::ostream &str) {
  
  str<<" List of Particle Data \n";
  str<<"      IDName";
  str<<std::setw(10)<<"kfc";
  str<<std::setw(14)<<"MASS[<kfc>]";
  str<<std::setw(16)<<"WIDTH[<kfc>]";
  str<<std::setw(16)<<"STABLE[<kfc>]";
  str<<std::setw(16)<<"MASSIVE[<kfc>]";
  str<<std::setw(16)<<"ACTIVE[<kfc>]\n";

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();
  
  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    kf_code fd(flav.Kfcode());
    // suppress pseudoparticle output
    while (fd/10) fd/=10;
    if (((!flav.IsHadron() && fd!=9) ||
	 (flav.Kfcode()>9900000 && flav.Kfcode()<9900099)) && 
	flav.Size()==1 && flav.Kfcode()!=0) {
      str<<std::setw(12)<<flav.IDName();
      str<<std::setw(10)<<flav.Kfcode();
      str<<std::setw(14)<<flav.Mass(true);
      str<<std::setw(16)<<flav.Width();
      str<<std::setw(16)<<flav.Stable();
      str<<std::setw(16)<<flav.IsMassive();
      str<<std::setw(15)<<flav.IsOn();
      str<<"\n";    
    }
  }
  str<<"\n";
}

void ATOOLS::OutputContainers(std::ostream &str) {
  
  str<<" List of Particle Containers \n";
  str<<"    IDName";
  str<<std::setw(8)<<"kfc";
  str<<std::setw(18)<<"Constituents\n";

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();
  
  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (!flav.IsHadron() && flav.Size()>1 && flav.Kfcode()!=0) {
      str<<std::setw(10)<<flav.IDName();
      str<<std::setw(8)<<flav.Kfcode();
      str<<std::setw(6)<<"{";
      for (unsigned int i=0;i<flav.Size();i++) {
	if (i!=flav.Size()-1) str<<flav[i].IDName()<<",";
	if (i==flav.Size()-1) str<<flav[i].IDName();
      }
      str<<"}\n";
    }
  }
  str<<"\n";
}

Mass_Selector::~Mass_Selector()
{
}
