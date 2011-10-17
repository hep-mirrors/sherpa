#include "METOOLS/Main/Polarization_Index.H"

#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;

void Polarization_Index::Init(const std::vector<int> &spins)
{
  m_spins=spins;
  if (m_spins.empty()) THROW(fatal_error,"No spin information");
  m_n=1;
  for(size_t i=0;i<m_spins.size();++i) m_n*=m_spins[i];
}

size_t Polarization_Index::operator()(const std::vector<int> &spins) const
{
  if(spins.size()!=m_spins.size())
    THROW(fatal_error,"Invalid size of spin vector");
  size_t mult=1, num=0;
  for(size_t i=0;i<spins.size();++i) {
    if(spins[i]<0 || spins[i]>m_spins[i])
      THROW(fatal_error,"Invalid spin index");
    num+=mult*spins[i];
    mult*=m_spins[i];
  }
  return num;
}

std::vector<int> Polarization_Index::operator()(size_t number) const
{
  std::vector<int> spins(m_spins.size());
  for(size_t i=0;i<m_spins.size();i++) {
    spins[i]=number%m_spins[i];
    number=(number-spins[i])/m_spins[i];
  }
  return spins;
}
