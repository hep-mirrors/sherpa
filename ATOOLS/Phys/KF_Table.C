#include "ATOOLS/Phys/KF_Table.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

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

  str<<"Particle data:\n";
  str<<Frame_Header(tablewidth);
  str<<Frame_Line(MyStrStream{}
     <<std::setw(9)<<"Name"
     <<std::setw(9)<<"Kf-code"
     <<std::setw(14)<<"Mass"
     <<std::setw(14)<<"Width"
     <<std::setw(9)<<"Stable"
     <<std::setw(9)<<"Massive"
     <<std::setw(9)<<"Active"
     <<std::setw(14)<<"Yukawa", tablewidth);
  str<<Frame_Separator(tablewidth);

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (flav.IsDiQuark() || flav.IsHadron()) continue;
    if (flav.Size()==1 && flav.Kfcode()!=0 && !flav.IsDummy()) {
      str<<Frame_Line(MyStrStream{}
         <<std::setw(9)<<flav.IDName()
         <<std::setw(9)<<flav.Kfcode()
         <<std::setw(14)<<flav.Mass(true)
         <<std::setw(14)<<flav.Width()
         <<std::setw(9)<<flav.Stable()
         <<std::setw(9)<<flav.IsMassive()
         <<std::setw(9)<<flav.IsOn()
         <<std::setw(14)<<flav.Yuk(), tablewidth);
    }
  }
  str<<Frame_Footer(tablewidth);
}

void ATOOLS::OutputContainers(std::ostream &str) {

  static int tablewidth {91};
  // There can be a lot of constituents, so we break the lines after a number
  // of constituents to prevent the output from becoming very wide.
  static int constituents_per_row {14};
  str<<"Particle containers:\n";
  str<<Frame_Header(tablewidth);
  str<<Frame_Line(MyStrStream{}
     <<std::setw(9)<<"Name"
     <<std::setw(9)<<"Kf-code"
     <<"  Constituents", tablewidth);
  str<<Frame_Separator(tablewidth);

  KFCode_ParticleInfo_Map::const_iterator kfit = s_kftable.begin();

  for (;kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (!flav.IsHadron() && flav.IsGroup() && flav.Kfcode()!=0) {
      for (int j=0; j < (flav.Size() - 1) / constituents_per_row + 1; j++) {
        MyStrStream line;
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
        str<<Frame_Line(line, tablewidth);
      }
    }
  }
  str<<Frame_Footer(tablewidth);
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
