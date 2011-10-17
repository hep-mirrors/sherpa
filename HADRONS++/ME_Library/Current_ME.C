#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Current_ME::Current_ME(const ATOOLS::Flavour_Vector& flavs,
                       const std::vector<int>& decayindices,
                       const std::string& name):
  HD_ME_Base(flavs,decayindices,name), p_c1(NULL), p_c2(NULL)
{
}

Current_ME::~Current_ME() {
  if (p_c1) delete p_c1;
  if (p_c2) delete p_c2;
}

void Current_ME::SetCurrent1(Current_Base* c1)
{
  p_c1=c1;
}

void Current_ME::SetCurrent2(Current_Base* c2)
{
  p_c2=c2;
}

void Current_ME::Calculate(const Vec4D_Vector& p, bool anti)
{
  p_c1->Calc(p, anti);
  p_c2->Calc(p, anti);
  
  std::vector<int> spins,spins1,spins2;
  for(size_t i=0;i<size();i++) {
    spins=GetSpinCombination(i);
    spins1.clear(); spins2.clear();
    for(size_t j=0;j<p_c1->DecayIndices().size();j++)
      spins1.push_back(spins[p_c1->DecayIndices()[j]]);
    for(size_t j=0;j<p_c2->DecayIndices().size();j++)
      spins2.push_back(spins[p_c2->DecayIndices()[j]]);
    // now we know the spin combinations in both currents
    // let's fill the results:
    (*this)[i]=m_factor*p_c1->Get(spins1)*p_c2->Get(spins2);
  }
}
