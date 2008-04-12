#include "Flavour.H"

#include "MathTools.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Data_Reader.H"
#include "Run_Parameter.H"
#include "Run_Parameter.H"
#include "Random.H"

namespace ATOOLS 
{
  KF_Table s_kftable;
}

using namespace ATOOLS;

Particle_Info::Particle_Info(const Particle_Info &info):
  m_kfc(info.m_kfc), m_mass(info.m_mass), m_yuk(info.m_yuk), 
  m_width(info.m_width),
  m_dg(info.m_dg), m_dm(info.m_dm), m_qoverp2(info.m_qoverp2), 
  m_icharge(info.m_icharge), m_isoweak(info.m_isoweak), 
  m_strong(info.m_strong), m_spin(info.m_spin), m_stable(info.m_stable), 
  m_masssign(info.m_masssign), m_dummy(info.m_dummy), m_majorana(info.m_majorana), 
  m_on(info.m_on), m_massive(info.m_massive), m_hadron(info.m_hadron),
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
 const bool stable,bool massive,const std::string &idname,
 const std::string &texname,const bool dummy):
  m_kfc(kfc), m_mass(mass), m_yuk(mass), m_width(width),
  m_dg(0.0), m_dm(0.0), m_qoverp2(1.0),
  m_icharge(icharge), m_isoweak(isoweak), m_strong(strong), m_spin(spin), 
  m_stable(stable), m_masssign(1), m_dummy(dummy), m_majorana(majorana), 
  m_on(on), m_massive(massive), m_hadron(0), m_idname(idname), 
  m_texname(texname)
{
  m_content.push_back(new Flavour(m_kfc));
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass,const double &width,
 const int icharge,const int isoweak,const int spin,const bool on,
 const bool stable,const std::string &idname,const std::string &texname):
  m_kfc(kfc), m_mass(mass), m_yuk(0.0), m_width(width),
  m_icharge(icharge), m_isoweak(isoweak), m_strong(0), m_spin(spin), 
  m_stable(stable), m_masssign(1), m_dummy(0), m_majorana(0), m_on(on), 
  m_massive(1), m_hadron(1), m_idname(idname), m_texname(texname) 
{
  m_content.push_back(new Flavour(m_kfc));
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
  m_kfc=(kf_code)abs(code-6);
  if (code==6) m_kfc=kf_gluon;
}

int Flavour::HepEvt() 
{
  switch (m_kfc) {
  case kf_a_0_1450:      return 10111;
  case kf_a_0_1450_plus: return m_anti?-10211:10211;
  case kf_f_0_1370:      return 10221;
  case kf_f_0_1710:      return 10331;
  case kf_a_0_980:       return 9000111;
  case kf_a_0_980_plus:  return m_anti?-9000211:9000211;
  case kf_f_0_980:       return 9010221;
  }
  return (long int)*this;
}

void Flavour::FromHepEvt(long int code) 
{
  m_anti=code<0;
  m_kfc=code=(kf_code)abs(code);
  switch (code) {
  case 10111:   m_kfc=kf_a_0_1450; break;
  case 10211:   m_kfc=kf_a_0_1450_plus; break;
  case 10221:   m_kfc=kf_f_0_1370; break;
  case 10331:   m_kfc=kf_f_0_1710; break;
  case 9000111: m_kfc=kf_a_0_980; break;
  case 9000211: m_kfc=kf_a_0_980_plus; break;
  case 9010221: m_kfc=kf_f_0_980; break; 
  case 91:      m_kfc=kf_cluster; break;
  case 92:      m_kfc=kf_string; break;
  }
}

std::string Flavour::TexName() const 
{
  if (SelfAnti()) return s_kftable[m_kfc]->m_texname;
  std::string name;
  if (m_anti) name="\\overline{";
  switch (m_kfc) {
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
      name+=s_kftable[m_kfc]->m_texname;
    }
    break;
  }
  if (m_anti) {
    name+="}";
    switch (m_kfc) {
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
  std::string name(s_kftable[m_kfc]->m_idname);
  if (m_kfc==kf_e || m_kfc==kf_mu || m_kfc==kf_tau) {
    name.erase(name.length()-1,1);
    if (IsAnti()) name+="+";
    else name+="-";      
  }
  else {
    if (m_kfc==kf_Hplus || m_kfc==kf_Wplus ||
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
  if(abs(m_kfc)>=1103&&abs(m_kfc)<=5505) {
    double help=abs(m_kfc)/100.0-int(abs(m_kfc)/100.0); 
    if(help<0.031) return true;
  }
  return false;
}

bool Flavour::IsBaryon() const 
{
  if (abs(m_kfc)%10000<1000) return false;
  return !IsDiQuark();
}

bool Flavour::IsB_Hadron() const 
{
  if (abs(m_kfc)<100)                            return 0;
  if (m_kfc-100*int(m_kfc/100)<10)                 return 0;
  if (abs((m_kfc-100*int(m_kfc/100))/10)==5)       return 1;
  if (abs((m_kfc-1000*int(m_kfc/1000))/100)==5)    return 1;
  if (abs((m_kfc-10000*int(m_kfc/10000))/1000)==5) return 1;
  return 0;
}

bool Flavour::IsC_Hadron() const 
{
  if (abs(m_kfc)<100)                            return 0;
  if (m_kfc-100*int(m_kfc/100)<10)                 return 0;
  if (abs((m_kfc-100*int(m_kfc/100))/10)==4)       return 1;
  if (abs((m_kfc-1000*int(m_kfc/1000))/100)==4)    return 1;
  if (abs((m_kfc-10000*int(m_kfc/10000))/1000)==4) return 1;
  return 0;
}

double Flavour::DiceLifeTime() const
{
  double proper_time = rpa.hBar() / Width();
  return -proper_time*log(1.-ran.Get());
}

std::ostream &ATOOLS::operator<<(std::ostream &os,const Flavour &fl)
{
  return os<<fl.IDName();
}

void ATOOLS::ParticleInit(const std::string &path)
{ 
  static bool initialized(false);
  if (initialized) return;
  initialized=true;
  kf_code kfc;
  int    charge,icharge,spin,strong,Majorana;
  bool   Take,stable,massive;
  double mass,width;
  std::string idname, texname, buffer;
  
  Data_Reader read(" ",";","!","=");
  read.AddWordSeparator("\t");
  read.SetAddCommandLine(false);
  read.SetInputPath(path);

  /*
  read.SetInputFile(rpa.gen.Variable("PARTICLE_DATA_FILE"));

  if (!read.OpenInFile()) {
    msg_Error()<<METHOD<<"(): File '"
	       <<rpa.gen.Variable("PARTICLE_DATA_FILE")
	       <<"' not found."<<std::endl;
    return;
  }
  msg_LogFile()<<"\n! +--------------------------------------+"<<std::endl;
  msg_LogFile()<<"! |               kf table               |"<<std::endl;
  msg_LogFile()<<"! +--------------------------------------+\n"<<std::endl;
  msg_LogFile()<<"! "<<buffer<<std::endl;
  */

  std::map<int,double> cdm, cdw;
  std::map<int,int> cia, cis, cim;
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.AddIgnore("[");
  dr.AddIgnore("]");
  std::vector<std::vector<double> > helpdvv;
  if (dr.MatrixFromFile(helpdvv,"MASS"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdm[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"WIDTH"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdw[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"ACTIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cia[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"STABLE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cis[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"MASSIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cim[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  
  std::vector<std::vector<std::string> > helpsvv;
  /*
  read.MatrixFromFile(helpsvv);
  for(size_t i(1);i<helpsvv.size();++i) {
    if (helpsvv[i].size()!=13) {
      msg_Error()<<METHOD<<"(): Inconsistent entry in line "<<i
  		 <<" of '"<<read.InputFile()<<"'."<<std::endl;
      continue;
   }
    kfc=ToType<int>(helpsvv[i][0]); 
    mass=ToType<double>(helpsvv[i][1]); width=ToType<double>(helpsvv[i][2]);
    charge=ToType<int>(helpsvv[i][3]); icharge=ToType<int>(helpsvv[i][4]);
    strong=ToType<int>(helpsvv[i][5]); spin=ToType<int>(helpsvv[i][6]);
    Majorana=ToType<int>(helpsvv[i][7]); Take=ToType<int>(helpsvv[i][8]);
    stable=ToType<int>(helpsvv[i][9]); massive=ToType<int>(helpsvv[i][10]);
    idname=helpsvv[i][11]; texname=helpsvv[i][12];
      std::map<int,double>::const_iterator dit;
      if ((dit=cdm.find(kfc))!=cdm.end()) mass=dit->second;
      if ((dit=cdw.find(kfc))!=cdw.end()) width=dit->second;
      std::map<int,int>::const_iterator iit;
      if ((iit=cia.find(kfc))!=cia.end()) Take=iit->second;
      if ((iit=cis.find(kfc))!=cis.end()) stable=iit->second;
      if ((iit=cim.find(kfc))!=cim.end()) massive=iit->second;
      msg_LogFile()<<"! "<<kfc<<" \t"<<mass<<" \t"<<width<<" \t"
  		   <<charge<<" \t"<<icharge<<" \t"<<strong<<" \t"<<spin
  		   <<" \t"<<Majorana<<" \t"<<Take<<" \t"<<stable<<" \t"
  		   <<massive<<" \t"<<idname<<"\t "<<texname<<std::endl;
     s_kftable[kfc] = new
  	Particle_Info(kfc,mass,width,charge,icharge,strong,spin,
  		      Majorana,Take,stable,massive,idname,texname);
  }
  msg_LogFile()<<std::endl;
  */
  read.SetFileBegin(std::vector<std::string>());
  read.SetFileEnd(std::vector<std::string>());
  read.SetInputFile(rpa.gen.Variable("HADRON_DATA_FILE"));

  read.RescanInFile();
  if (!read.OpenInFile()) {
    msg_Error()<<METHOD<<"(): File '"
	       <<rpa.gen.Variable("HADRON_DATA_FILE")
	       <<"' not found."<<std::endl;
  }
  else {
    msg_LogFile()<<"! "<<buffer<<std::endl;
    read.MatrixFromFile(helpsvv);
    for(size_t i(1);i<helpsvv.size();++i) {
      if (helpsvv[i].size()!=9) {
	msg_Error()<<METHOD<<"(): Inconsistent entry in line "<<i
		   <<" of '"<<read.InputFile()<<"'."<<std::endl;
	continue;
      }
      kfc=ToType<int>(helpsvv[i][0]); 
      mass=ToType<double>(helpsvv[i][1]); width=ToType<double>(helpsvv[i][2]);
      charge=ToType<int>(helpsvv[i][3]); icharge=ToType<int>(helpsvv[i][4]);
      spin=ToType<int>(helpsvv[i][5]); Take=ToType<int>(helpsvv[i][6]);
      stable=ToType<int>(helpsvv[i][7]); idname=helpsvv[i][8];
      texname=idname;
	std::map<int,double>::const_iterator dit;
	if ((dit=cdm.find(kfc))!=cdm.end()) mass=dit->second;
	if ((dit=cdw.find(kfc))!=cdw.end()) width=dit->second;
	std::map<int,int>::const_iterator iit;
	if ((iit=cia.find(kfc))!=cia.end()) Take=iit->second;
	if ((iit=cis.find(kfc))!=cis.end()) stable=iit->second;
	msg_LogFile()<<"! "<<kfc<<" \t"<<mass<<" \t"<<width<<" \t"
		     <<charge<<" \t"<<icharge<<" \t"<<spin<<" \t"
		     <<Take<<" \t"<<stable<<" \t"<<idname<<"\t "
		     <<texname<<std::endl;
	s_kftable[kfc] = new
	  Particle_Info(kfc,mass,width,charge,icharge, 
			spin,Take,stable,idname,texname);
        if(kfc==kf_pi||kfc==kf_eta||kfc==kf_eta_prime_958||kfc==kf_eta_c_1S||
           kfc==kf_eta_b||kfc==kf_rho_770||kfc==kf_omega_782||
           kfc==kf_phi_1020||kfc==kf_J_psi_1S||kfc==kf_Upsilon_1S||
           kfc==kf_a_2_1320||kfc==kf_f_2_1270||kfc==kf_f_2_prime_1525||
           kfc==kf_chi_c2_1P||kfc==kf_chi_b2_1P) {
          s_kftable[kfc]->m_majorana=-1;
        }
    }
    msg_LogFile()<<std::endl;
  }

  s_kftable[kf_none] = new
    Particle_Info(kf_none,-1,0,0,0,0,0,-1,0,1,0,"no_particle","no particle");
  s_kftable[kf_cluster] = new
    Particle_Info(kf_cluster,0.,0.,0,0,0,0,0,1,1,0,"cluster","cluster");
  s_kftable[kf_string] = new
    Particle_Info (kf_string,0.,0.,0,0,0,0,0,1,1,0,"string","string");
  s_kftable[kf_seaquark] = new
    Particle_Info(kf_seaquark,0.,0.,0,0,1,1,0,1,1,0,"seaquark","seaquark");
  s_kftable[kf_bjet] = new
    Particle_Info(kf_bjet,0.,0.,0,0,1,2,0,1,1,0,"bj","bjet");

  s_kftable[kf_fermion] = new
    Particle_Info(kf_fermion,0.,0., 0,0,0,1,0,1,1,0,"fermion","fermion");
  s_kftable[kf_jet] = new
    Particle_Info(kf_jet,0.,0.,0,0,1, 2,0,1,1,0,"j","jet");
  s_kftable[kf_quark] = new
    Particle_Info(kf_quark,0.,0.,0, 0,1,1,0,1,1,0,"Q","Quark");
  s_kftable[kf_lepton] = new
    Particle_Info(kf_lepton,0.,0.,-3,-1,0,1,0,1,1,0,"lepton","lepton");
  s_kftable[kf_neutrino] = new
    Particle_Info(kf_neutrino,0.,0.,0,1,0, 1,0,1,1,0,"neutrino","neutrino");
  s_kftable[kf_fermion]->Clear();
  s_kftable[kf_jet]->Clear();
  s_kftable[kf_quark]->Clear();
  s_kftable[kf_lepton]->Clear();
  s_kftable[kf_neutrino]->Clear();
  for (int i=1;i<7;i++) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0) {
      s_kftable[kf_jet]->Add(addit);
      s_kftable[kf_jet]->Add(addit.Bar());
      s_kftable[kf_quark]->Add(addit);
      s_kftable[kf_quark]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
  s_kftable[kf_jet]->Add(Flavour(kf_gluon));
  for (int i=11;i<17;i+=2) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0) {
      s_kftable[kf_lepton]->Add(addit);
      s_kftable[kf_lepton]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
  for (int i=12;i<18;i+=2) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0) {
      s_kftable[kf_neutrino]->Add(addit);
      s_kftable[kf_neutrino]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
}
  
