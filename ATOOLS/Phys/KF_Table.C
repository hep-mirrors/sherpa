#include "ATOOLS/Phys/KF_Table.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <iomanip>

namespace ATOOLS
{
  KF_Table s_kftable;
}

using namespace ATOOLS;

void ATOOLS::OutputHadrons(std::ostream &str) {

  str<<"List of Hadron data \n";
  str<<std::setw(16)<<"IDName";
  str<<std::setw(8)<<"kfc";
  str<<std::setw(16)<<"Mass";
  str<<std::setw(16)<<"Width";
  str<<std::setw(9)<<"Stable";
  str<<std::setw(9)<<"Active";
  str<<'\n';

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if ((flav.IsHadron() || flav.IsDiQuark())
        && flav.Size()==1 && flav.Kfcode()!=0) {
      str<<std::setw(16)<<flav.IDName();
      str<<std::setw(8)<<flav.Kfcode();
      str<<std::setw(16)<<flav.HadMass();
      str<<std::setw(16)<<flav.Width();
      str<<std::setw(9)<<flav.Stable();
      str<<std::setw(9)<<flav.IsOn();
      str<<"\n";
    }
  }
}

void ATOOLS::OutputParticles(std::ostream &str) {

  static int tablewidth {91};
  MyStrStream line;

  str<<"Particle data:\n";
  str<<Frame_Header(tablewidth);
  line<<std::setw(9)<<"Name"
      <<std::setw(9)<<"Kf-code"
      <<std::setw(14)<<"Mass"
      <<std::setw(14)<<"Width"
      <<std::setw(9)<<"Stable"
      <<std::setw(9)<<"Massive"
      <<std::setw(9)<<"Active"
      <<std::setw(14)<<"Yukawa";
  str<<Frame_Line(line.str(), tablewidth);
  str<<Frame_Separator(tablewidth);

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (flav.IsDiQuark() || flav.IsHadron()) continue;
    if (flav.Size()==1 && flav.Kfcode()!=0 && !flav.IsDummy()) {
      line.str("");
      line<<std::setw(9)<<flav.IDName()
          <<std::setw(9)<<flav.Kfcode()
          <<std::setw(14)<<flav.Mass(true)
          <<std::setw(14)<<flav.Width()
          <<std::setw(9)<<flav.Stable()
          <<std::setw(9)<<flav.IsMassive()
          <<std::setw(9)<<flav.IsOn()
          <<std::setw(14)<<flav.Yuk();
      str<<Frame_Line(line.str(), tablewidth);
    }
  }
  str<<Frame_Footer(tablewidth);
}

void ATOOLS::OutputContainers(std::ostream &str) {

  static int tablewidth {91};
  MyStrStream line;
  // There can be a lot of constituents, so we break the lines after a number
  // of constituents to prevent the output from becoming very wide.
  static int constituents_per_row {14};
  str<<"Particle containers:\n";
  str<<Frame_Header(tablewidth);
  line<<std::setw(9)<<"Name"
      <<std::setw(9)<<"Kf-code"
      <<"  Constituents";
  str<<Frame_Line(line.str(), tablewidth);
  str<<Frame_Separator(tablewidth);

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (!flav.IsHadron() && flav.IsGroup() && flav.Kfcode()!=0 && flav.Size()!=0) {
      for (int j=0; j < (flav.Size() - 1) / constituents_per_row + 1; j++) {
        line.str("");
        if (j==0) {
          line<<std::setw(9)<<flav.IDName();
          line<<std::setw(9)<<flav.Kfcode();
        } else {
          line<<std::setw(9)<<"";
          line<<std::setw(9)<<"";
        }
        line<<"  ";
        for (unsigned int i=j*constituents_per_row;i<Min((int)flav.Size(),(j+1)*constituents_per_row);i++) {
          if (i!=flav.Size()-1) line<<flav[i].IDName()<<", ";
          if (i==flav.Size()-1) line<<flav[i].IDName();
        }
        str<<Frame_Line(line.str(), tablewidth);
      }
    }
  }
  str<<Frame_Footer(tablewidth);
}

void ATOOLS::AddParticle(kf_code kfc, double mass, double radius, double width,
                         int icharge, int strong, int spin, int majorana,
                         bool on, int stable, bool massive,
                         const std::string &idname, const std::string &antiname,
                         const std::string &texname,
                         const std::string &antitexname,
                         bool dummy, bool group)
{
  if (s_kftable.find(kfc) != s_kftable.end()) {
    THROW(fatal_error, "Particle already added.");
  }
  s_kftable[kfc] = new Particle_Info(
      kfc, mass, radius, width, icharge, strong, spin, majorana, on, stable,
      massive, idname, antiname, texname, antitexname, dummy, group);
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  const std::string key {ToString(kfc)};
  s_kftable[kfc]->m_mass = s_kftable[kfc]->m_hmass =
      pdata[key]["Mass"].SetDefault(mass).Get<double>();
  s_kftable[kfc]->m_width =
      pdata[key]["Width"].SetDefault(width).Get<double>();
  s_kftable[kfc]->m_on =
      pdata[key]["Active"].SetDefault(on).Get<bool>();
  s_kftable[kfc]->m_stable =
      pdata[key]["Stable"].SetDefault(stable).Get<int>();
  s_kftable[kfc]->m_massive =
      pdata[ToString(kfc)]["Massive"].SetDefault(massive).Get<bool>();
  s_kftable[kfc]->m_icharge =
      pdata[key]["IntCharge"].SetDefault(icharge).Get<int>();
  s_kftable[kfc]->m_strong =
      pdata[key]["StrongCharge"].SetDefault(strong).Get<int>();
  s_kftable[kfc]->m_yuk =
      pdata[key]["Yukawa"].SetDefault(s_kftable[kfc]->m_yuk).Get<double>();
  s_kftable[kfc]->m_priority =
      pdata[key]["Priority"].SetDefault(s_kftable[kfc]->m_priority).Get<int>();
}

void ATOOLS::AddOrUpdateParticle(kf_code kfc, double mass, double radius,
                                 double width, int icharge, int strong,
                                 int spin, int majorana, bool on, int stable,
                                 bool massive, const std::string &idname,
                                 const std::string &antiname,
                                 const std::string &texname,
                                 const std::string &antitexname, bool dummy,
                                 bool group)
{
  if (s_kftable.find(kfc) == s_kftable.end()) {
    s_kftable[kfc] = new Particle_Info(
        kfc, mass, radius, width, icharge, strong, spin, majorana, on, stable,
        massive, idname, antiname, texname, antitexname, dummy, group);
  } else {
    s_kftable[kfc]->m_radius = radius;
    s_kftable[kfc]->m_spin = spin;
    s_kftable[kfc]->m_majorana = majorana;
    s_kftable[kfc]->m_dummy = dummy;
    s_kftable[kfc]->m_isgroup = group;
    s_kftable[kfc]->m_idname = idname;
    s_kftable[kfc]->m_antiname = antiname;
    s_kftable[kfc]->m_texname = texname;
    s_kftable[kfc]->m_antitexname = antitexname;
  }
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  const std::string key {ToString(kfc)};
  s_kftable[kfc]->m_mass = s_kftable[kfc]->m_hmass =
      pdata[key]["Mass"].GetScalarWithOtherDefault(mass);
  s_kftable[kfc]->m_width =
      pdata[key]["Width"].GetScalarWithOtherDefault(width);
  s_kftable[kfc]->m_on =
      pdata[key]["Active"].GetScalarWithOtherDefault(on);
  s_kftable[kfc]->m_stable =
      pdata[key]["Stable"].GetScalarWithOtherDefault(stable);
  s_kftable[kfc]->m_massive =
      pdata[key]["Massive"].GetScalarWithOtherDefault(massive);
  s_kftable[kfc]->m_icharge =
      pdata[key]["IntCharge"].GetScalarWithOtherDefault(icharge);
  s_kftable[kfc]->m_strong =
      pdata[key]["StrongCharge"].GetScalarWithOtherDefault(strong);
}

void ATOOLS::AddParticle(kf_code kfc, double mass, double radius, double width,
                         int icharge, int spin, int majorana, bool on,
                         int stable, const std::string &idname,
                         const std::string &texname)
{
  if (s_kftable.find(kfc) != s_kftable.end()) {
    THROW(fatal_error, "Particle already added.");
  }
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  const std::string key {ToString(kfc)};
  s_kftable[kfc] = new Particle_Info(kfc, mass, radius, width, icharge, spin,
                                     on, stable, idname, texname);
  s_kftable[kfc]->m_mass = s_kftable[kfc]->m_hmass =
      pdata[key]["Mass"].SetDefault(mass).Get<double>();
  s_kftable[kfc]->m_width =
      pdata[key]["Width"].SetDefault(width).Get<double>();
  s_kftable[kfc]->m_on =
      pdata[key]["Active"].SetDefault(on).Get<bool>();
  s_kftable[kfc]->m_stable =
      pdata[key]["Stable"].SetDefault(stable).Get<int>();
  s_kftable[kfc]->m_massive =
      pdata[key]["Massive"].SetDefault(s_kftable[kfc]->m_massive).Get<bool>();
  s_kftable[kfc]->m_icharge =
      pdata[key]["IntCharge"].SetDefault(icharge).Get<int>();
  s_kftable[kfc]->m_strong =
      pdata[key]["StrongCharge"].SetDefault(s_kftable[kfc]->m_strong).Get<int>();
  s_kftable[kfc]->m_yuk =
      pdata[key]["Yukawa"].SetDefault(s_kftable[kfc]->m_yuk).Get<double>();
  s_kftable[kfc]->m_priority =
      pdata[key]["Priority"].SetDefault(s_kftable[kfc]->m_priority).Get<int>();
}

void ATOOLS::AddParticle(kf_code kfc, double mass, double radius, double width,
                         int icharge, int spin, bool on, int stable,
                         const std::string &idname,
                         const std::string &texname)
{
  if (s_kftable.find(kfc) != s_kftable.end()) {
    THROW(fatal_error, "Particle already added.");
  }
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  const std::string key {ToString(kfc)};
  s_kftable[kfc]=new Particle_Info(kfc, mass, radius, width, icharge, spin, on, stable,
                                   idname, texname);
  s_kftable[kfc]->m_mass = s_kftable[kfc]->m_hmass =
      pdata[key]["Mass"].SetDefault(mass).Get<double>();
  s_kftable[kfc]->m_width =
      pdata[key]["Width"].SetDefault(width).Get<double>();
  s_kftable[kfc]->m_on =
      pdata[key]["Active"].SetDefault(on).Get<bool>();
  s_kftable[kfc]->m_stable =
      pdata[key]["Stable"].SetDefault(stable).Get<int>();
  s_kftable[kfc]->m_massive =
      pdata[key]["Massive"].SetDefault(s_kftable[kfc]->m_massive).Get<bool>();
  s_kftable[kfc]->m_icharge =
      pdata[key]["IntCharge"].SetDefault(icharge).Get<int>();
  s_kftable[kfc]->m_strong =
      pdata[key]["StrongCharge"].SetDefault(s_kftable[kfc]->m_strong).Get<int>();
  s_kftable[kfc]->m_yuk =
      pdata[key]["Yukawa"].SetDefault(s_kftable[kfc]->m_yuk).Get<double>();
  s_kftable[kfc]->m_priority =
      pdata[key]["Priority"].SetDefault(s_kftable[kfc]->m_priority).Get<int>();
}

void ATOOLS::AddParticle(kf_code kfc, double mass, double radius, int icharge,
                         int spin, int formfactor, const std::string &idname,
                         const std::string &texname)
{
  if (s_kftable.find(kfc) != s_kftable.end()) {
    THROW(fatal_error, "Particle already added.");
  }
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  const std::string key {ToString(kfc)};
  s_kftable[kfc] = new Particle_Info(kfc, mass, radius, icharge, spin,
                                     formfactor, idname, texname);
  s_kftable[kfc]->m_mass = s_kftable[kfc]->m_hmass =
      pdata[key]["Mass"].SetDefault(mass).Get<double>();
  s_kftable[kfc]->m_width =
      pdata[key]["Width"].SetDefault(s_kftable[kfc]->m_width).Get<double>();
  s_kftable[kfc]->m_on =
      pdata[key]["Active"].SetDefault(s_kftable[kfc]->m_on).Get<bool>();
  s_kftable[kfc]->m_stable =
      pdata[key]["Stable"].SetDefault(s_kftable[kfc]->m_stable).Get<int>();
  s_kftable[kfc]->m_massive =
      pdata[key]["Massive"].SetDefault(s_kftable[kfc]->m_massive).Get<bool>();
  s_kftable[kfc]->m_icharge =
      pdata[key]["IntCharge"].SetDefault(icharge).Get<int>();
  s_kftable[kfc]->m_strong =
      pdata[key]["StrongCharge"].SetDefault(s_kftable[kfc]->m_strong).Get<int>();
  s_kftable[kfc]->m_yuk =
      pdata[key]["Yukawa"].SetDefault(s_kftable[kfc]->m_yuk).Get<double>();
  s_kftable[kfc]->m_priority =
      pdata[key]["Priority"].SetDefault(s_kftable[kfc]->m_priority).Get<int>();
}

void ATOOLS::ClearParticles()
{
  for (KF_Table::const_iterator kfit(s_kftable.begin());kfit!=s_kftable.end();++kfit)
    delete kfit->second;
  s_kftable.clear();
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
