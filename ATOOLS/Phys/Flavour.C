#include "ATOOLS/Phys/Flavour.H"

#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

Particle_Info::Particle_Info(const Particle_Info &info):
  m_kfc(info.m_kfc), m_mass(info.m_mass), m_hmass(info.m_hmass), m_radius(info.m_radius),
  m_yuk(info.m_yuk), m_width(info.m_width),
  m_dg(info.m_dg), m_dm(info.m_dm), m_qoverp2(info.m_qoverp2), 
  m_icharge(info.m_icharge),
  m_strong(info.m_strong), m_resummed(info.m_resummed), m_priority(info.m_priority),
  m_spin(info.m_spin), m_stable(info.m_stable), 
  m_masssign(info.m_masssign), m_dummy(info.m_dummy), m_majorana(info.m_majorana), 
  m_formfactor(0), m_on(info.m_on), m_massive(info.m_massive), m_hadron(info.m_hadron),
  m_isgroup(info.m_isgroup), m_idname(info.m_idname), m_antiname(info.m_antiname),
  m_texname(info.m_texname), m_antitexname(info.m_antitexname)
{
  m_content.resize(info.m_content.size());
  for (size_t i(0);i<info.m_content.size();++i) 
    m_content[i] = new Flavour(*info.m_content[i]); 
}

Particle_Info::Particle_Info
(const kf_code &kfc, const double &mass, const double &radius, const double &width,
 const int icharge, const int strong,
 const int spin, const int majorana, const bool on,
 const int stable, bool massive, const std::string &idname,
 const std::string &antiname, const std::string& texname,
 const std::string &antitexname, const bool dummy, const bool isgroup):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_radius(radius), m_yuk(-1.0), m_width(width),
  m_dg(0.0), m_dm(0.0), m_qoverp2(1.0), m_icharge(icharge),
  m_strong(strong), m_resummed(0), m_priority(0), m_spin(spin), 
  m_stable(stable), m_masssign(1), m_dummy(dummy), m_majorana(majorana), 
  m_formfactor(0), m_on(on), m_massive(massive), m_hadron(0), 
  m_isgroup(isgroup), m_idname(idname), m_antiname(antiname),
  m_texname(texname), m_antitexname(antitexname)  
{
  m_content.push_back(new Flavour(*this));
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass,const double &radius,const double &width,
 const int icharge,const int spin,const bool on,
 const int stable,const std::string &idname,const std::string &texname):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_radius(radius), m_yuk(-1.0), m_width(width),
  m_dg(0.0), m_dm(0.0), m_qoverp2(1.0),
  m_icharge(icharge), m_strong(0), m_resummed(0), m_priority(0), 
  m_spin(spin), m_stable(stable), m_masssign(1), m_dummy(0), m_majorana(0), 
  m_formfactor(0), m_on(on), m_massive(1), m_hadron(1), m_isgroup(0), 
  m_idname(idname), m_texname(texname)
{
  m_antiname=m_idname+"b";
  m_antitexname="\\overline{"+m_antiname+"}";
  m_content.push_back(new Flavour(*this));
}

Particle_Info::Particle_Info
(const kf_code &kfc,const double &mass, const double &radius,const int icharge, const int spin,
 const int formfactor, const std::string &idname, const std::string &antiname):
  m_kfc(kfc), m_mass(mass), m_hmass(mass), m_radius(radius), m_yuk(-1.0), m_width(0),
  m_dg(0.0), m_dm(0.0), m_qoverp2(1.0),
  m_icharge(icharge), m_strong(0), m_resummed(0), m_priority(0), m_spin(spin), 
  m_stable(1), m_masssign(1), m_dummy(0), m_majorana(0), 
  m_formfactor(formfactor), m_on(1), m_massive(1), m_hadron(1), m_isgroup(0), 
  m_idname(idname), m_antiname(antiname)
{
  m_antiname=m_idname+"b";
  m_antitexname="\\overline{"+m_antiname+"}";
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
  }
}

void Particle_Info::Add(const Flavour &fl)
{ 
  if (m_mass>=0.0) {
    if (m_content.empty()) {
      m_mass=fl.Mass();
      m_massive=fl.IsMassive();
    }
    else {
      if (m_mass!=fl.Mass()) {
	msg_Error()<<METHOD<<"(): m_{"<<m_idname<<"} = "<<m_mass
		   <<" vs. m_{"<<fl<<"} = "<<fl.Mass(true)<<std::endl;
	THROW(critical_error,"Inconsistent input");
      }
      if (m_massive!=fl.IsMassive()) {
	msg_Error()<<METHOD<<"(): massive_{"<<m_idname<<"} = "<<m_massive
		   <<" vs. massive_{"<<fl<<"} = "<<fl.IsMassive()<<std::endl;
	THROW(critical_error,"Inconsistent input");
      }
    }
  }
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

Flavour::Flavour(Particle_Info& info, bool anti):
p_info(&info),
m_anti(0)
{
  if (anti && p_info->m_majorana==0) m_anti=anti;
}

Flavour::Flavour(long int kfc):
p_info(NULL),
m_anti(0)
{
  KFCode_ParticleInfo_Map::iterator it(s_kftable.find(std::abs(kfc)));
  if (it!=s_kftable.end()) p_info=it->second; else return;
  if (kfc<0 && p_info->m_majorana==0) m_anti=1;
}

Flavour::Flavour(kf_code kfc, bool anti):
p_info(NULL),
m_anti(0)
{
  KFCode_ParticleInfo_Map::iterator it(s_kftable.find(kfc));
  if (it!=s_kftable.end()) p_info=it->second; else return;
  if (anti && p_info->m_majorana==0) m_anti=anti;
}

Flavour::Flavour(const Flavour& fl):
p_info(fl.p_info),
m_anti(fl.m_anti)
{
}

void Flavour::InitializeParticleInfo(kf_code kfc)
{
  const auto it = s_kftable.find(kfc);
  if (it != s_kftable.end())
    p_info = it->second;
  else
    THROW(fatal_error, "Unknown particle code " + ToString(kfc));
}

double Flavour::ISSymmetryFactor(const ATOOLS::Flavour_Vector& flavs)
{
  double sf(1.0);
  for (ATOOLS::Flavour_Vector::const_iterator 
	 it=flavs.begin(); it!=flavs.end(); ++it)
    {
      double pols(2.0*it->Spin()+1.0);
      if (it->IntSpin()==2 &&
	  it->Mass()==0.0) pols=2.0;
      sf*=pols;
      if (it->Strong())sf*=abs(it->StrongCharge());
    }
  return sf;
}

double Flavour::FSSymmetryFactor(const ATOOLS::Flavour_Vector& flavs)
{
  double sf(1.0);
  std::map<ATOOLS::Flavour,size_t> fc;
  for (ATOOLS::Flavour_Vector::const_iterator 
	 it=flavs.begin(); it!=flavs.end(); ++it)
    fc[*it] = 0;
  for (ATOOLS::Flavour_Vector::const_iterator 
	 it=flavs.begin(); it!=flavs.end(); ++it)
    fc[*it] +=1;
  for (std::map<ATOOLS::Flavour,size_t>::const_iterator 
	 it(fc.begin()); it!=fc.end(); ++it) 
    sf*=Factorial(it->second);
  return sf;
}

std::string Flavour::TexName() const 
{
  if (!IsHadron()) return m_anti?p_info->m_antitexname:p_info->m_texname;
  std::string name, idname(IDName());
  bool barit(false);
  if (m_anti && (!SelfAnti()) && IsHadron()) { 
    if (idname.find("++")!=std::string::npos ||
	idname.find("+")!=std::string::npos  ||
	idname.find("-")!=std::string::npos  ||
	idname.find("--")!=std::string::npos) barit = true;
    else name="\\bar ";
  }

  switch (Kfcode()) {
  case kf_pi : {name+=std::string("\\pi^{0}");break;}
  case kf_K : {name+=std::string("K^{0}");break;}
  case kf_K_L : {name+=std::string("K_{L}");break;}
  case kf_K_S : {name+=std::string("K_{S}");break;}
  case kf_D : {name+=std::string("D^{0}"); break;}
  case kf_B : {name+=std::string("B^{0}"); break;}
  default :
    name+=IDName();
    if (barit) {
      if (name.find("++")!=std::string::npos)      name=StringReplace(name, "++", "--");
      else if (name.find("+")!=std::string::npos)  name=StringReplace(name, "+", "-");
      else if (name.find("--")!=std::string::npos) name=StringReplace(name, "--", "++");
      else if (name.find("-")!=std::string::npos)  name=StringReplace(name, "-", "+");
    }
    
    if (name.find("++")!=std::string::npos)      name=StringReplace(name, "++", "^{++}");
    else if (name.find("--")!=std::string::npos) name=StringReplace(name, "--", "^{--}");
    else if (name.find("*+")!=std::string::npos) name=StringReplace(name, "*+", "^{*+}");
    else if (name.find("*-")!=std::string::npos) name=StringReplace(name, "*-", "^{*-}");
    else if (name.find("'+")!=std::string::npos) name=StringReplace(name, "'+", "^{\\prime +}");
    else if (name.find("'-")!=std::string::npos) name=StringReplace(name, "'-", "^{\\prime -}");
    else if (name.find("*")!=std::string::npos)  name=StringReplace(name, "*", "^{*}");
    else if (name.find("'")!=std::string::npos)  name=StringReplace(name, "'", "^{\\prime }");
    else if (name.find("+")!=std::string::npos)  name=StringReplace(name, "+", "^{+}");
    else if (name.find("-")!=std::string::npos)  name=StringReplace(name, "-", "^{-}");
    name=StringReplace(name, "_fict", "[\\mathrm{fict.}]");
    name=StringReplace(name, "chi", "\\chi ");
    name=StringReplace(name, "nu", "\\nu");
    name=StringReplace(name, "mu", "\\mu");
    name=StringReplace(name, "tau", "\\tau");
    name=StringReplace(name, "\\nu_eb", "\\bar\\nu_e");
    name=StringReplace(name, "\\nu_\\mub", "\\bar\\nu_\\mu");
    name=StringReplace(name, "\\nu_\\taub", "\\bar\\nu_\\tau");
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
    name=StringReplace(name, "U\\psi lon", "\\Upsilon ");
    name=StringReplace(name, " _", "_");
    if (IsAnti() && name.back()=='b') name.pop_back();
    break;
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
  while ((pos=name.find("~"))!=std::string::npos) name.replace(pos,1,"tilde");
  return name;
}

// Shell names used to recover decay data file
// names. Used by HADRONS::Hadron_Decay_Table
std::string Flavour::LegacyShellName() const
{
  if(!m_anti){
    switch (Kfcode()){
    case      kf_d:           return "d";
    case      kf_u:           return "u";
    case      kf_s:           return "s";
    case      kf_c:           return "c";
    case      kf_b:           return "b";
    case      kf_t:           return "t";
    case      kf_e:           return "e-";
    case      kf_nue:         return "nu_e";
    case      kf_mu:          return "mu-";
    case      kf_numu:        return "nu_mu";
    case      kf_tau:         return "tau-";
    case      kf_nutau:       return "nu_tau";
    case      kf_gluon:       return "G";
    case      kf_photon:      return "P";
    case      kf_Z:           return "Z";
    case      kf_Wplus:       return "W+";
    case      kf_h0:          return "h0";
    case      kf_phiplus:     return "phi+";
    case      kf_chi:         return "chi";
    }
  }
  else{
    switch (Kfcode()){
    case      kf_d:           return "db";
    case      kf_u:           return "ub";
    case      kf_s:           return "sb";
    case      kf_c:           return "cb";
    case      kf_b:           return "bb";
    case      kf_t:           return "tb";
    case      kf_e:           return "e+";
    case      kf_nue:         return "nu_eb";
    case      kf_mu:          return "mu+";
    case      kf_numu:        return "nu_mub";
    case      kf_tau:         return "tau+";
    case      kf_nutau:       return "nu_taub";
    case      kf_gluon:       return "G";
    case      kf_photon:      return "P";
    case      kf_Z:           return "Z";
    case      kf_Wplus:       return "W-";
    case      kf_h0:          return "h0";
    case      kf_phiplus:     return "phi-";
    case      kf_chi:         return "chi";
    }
  }
  return ShellName();
}

std::string Flavour::IDName() const 
{
  return m_anti?p_info->m_antiname:p_info->m_idname;
}

bool Flavour::IsDiQuark() const 
{
  if(Kfcode() >= 1103 && Kfcode() <= 5505) {
    double help = Kfcode()/100.0 - int(Kfcode()/100.0); 
    if(help<0.031) return true;
  }
  return false;
}

bool Flavour::IsBaryon() const 
{
  if (Kfcode() % 10000 < 1000) return false;
  return !IsDiQuark();
}

bool Flavour::IsMeson() const 
{
  if (Kfcode() % 1000 < 100) return false;
  return !IsDiQuark();
}

bool Flavour::IsNucleon() const
{
  return (Kfcode()==2212 || Kfcode()==2112);
}

bool Flavour::IsB_Hadron() const 
{
  if (Kfcode() < 100)                               return 0;
  if (Kfcode()-100*int(Kfcode()/100)<10)            return 0;
  if ((Kfcode()-100*int(Kfcode()/100))/10==5)       return 1;
  if ((Kfcode()-1000*int(Kfcode()/1000))/100==5)    return 1;
  if ((Kfcode()-10000*int(Kfcode()/10000))/1000==5) return 1;
  return 0;
}

bool Flavour::IsC_Hadron() const 
{
  if (Kfcode() < 100)                               return 0;
  if (Kfcode()-100*int(Kfcode()/100)<10)            return 0;
  if ((Kfcode()-100*int(Kfcode()/100))/10==4)       return 1;
  if ((Kfcode()-1000*int(Kfcode()/1000))/100==4)    return 1;
  if ((Kfcode()-10000*int(Kfcode()/10000))/1000==4) return 1;
  return 0;
}

double Flavour::GenerateLifeTime() const
{
  double proper_time = rpa->hBar() / Width();
  return -proper_time*log(1.-ran->Get());
}

double Flavour::RelBWMass(const double& min, const double& max,
                          double peak, double width) const
{
  if (peak<0.0) peak=Mass(true);
  if (width<0.0) width=Width();
  if( peak<1.e-6 || width/peak < 1.e-8) return peak;
  double random = ran->Get();
  double peak2 = peak*peak;
  double mw    = peak*width;
  double s;
  if (min==0.0 && max==std::numeric_limits<double>::max()) {
    s = peak2+mw*tan(M_PI*(random-0.5));
  }
  else {
    double smin = sqr(min); double smax = sqr(max);
    double ymax=atan((smin-peak2)/mw);
    double ymin=atan((smax-peak2)/mw);
    s = peak2+mw*tan(ymin + random*(ymax-ymin));
  }
  return sqrt(s);
}

bool Flavour::IsStable() const
{
  if (p_info->m_stable==0) return false;
  if (p_info->m_stable==1) return true;
  if (p_info->m_stable==2 && !IsAnti()) return true;
  if (p_info->m_stable==3 && IsAnti()) return true;
  return false;
}

Flavour Flavour::IsoWeakPartner() const
{
  if (IsoWeak() != 0) {
    auto code = Kfcode();
    if (code % 2 == 0)
      --code;
    else
      ++code;
    return Flavour(code, m_anti);
  }
  return *this;
}

Flavour Flavour::GoldstoneBosonPartner() const
{
  auto code = Kfcode();
  if (code == kf_Z)
    code = kf_chi;
  else if (code == kf_Wplus)
    code = kf_phiplus;
  return Flavour(code, m_anti);
}

std::ostream &ATOOLS::operator<<(std::ostream &os,const Flavour &fl)
{
  return os<<fl.IDName();
}

Mass_Selector::~Mass_Selector()
{
}
