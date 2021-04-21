#include "ATOOLS/Phys/KF_Table.H"

#include "ATOOLS/Org/Message.H"

#include <iomanip>

namespace ATOOLS
{
  KF_Table s_kftable;
}

using namespace ATOOLS;

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
  str<<std::setw(8)<<"IDName";
  str<<std::setw(8)<<"kfc";
  str<<std::setw(13)<<"MASS[<kfc>]";
  str<<std::setw(15)<<"WIDTH[<kfc>]";
  str<<std::setw(15)<<"STABLE[<kfc>]";
  str<<std::setw(15)<<"MASSIVE[<kfc>]";
  str<<std::setw(15)<<"ACTIVE[<kfc>]";
  str<<std::setw(16)<<"YUKAWA[<kfc>]\n";

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (flav.IsDiQuark() || flav.IsHadron()) continue;
    if (flav.Size()==1 && flav.Kfcode()!=0 && !flav.IsDummy()) {
      str<<std::setw(8)<<flav.IDName();
      str<<std::setw(8)<<flav.Kfcode();
      str<<std::setw(13)<<flav.Mass(true);
      str<<std::setw(15)<<flav.Width();
      str<<std::setw(15)<<flav.Stable();
      str<<std::setw(15)<<flav.IsMassive();
      str<<std::setw(15)<<flav.IsOn();
      str<<std::setw(15)<<flav.Yuk();
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
    if (!flav.IsHadron() && flav.IsGroup() && flav.Kfcode()!=0) {
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
