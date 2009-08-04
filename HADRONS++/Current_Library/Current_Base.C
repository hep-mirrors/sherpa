#include "HADRONS++/Current_Library/Current_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::Current_Base
#define PARAMETER_TYPE HADRONS::Flavour_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;
using namespace HADRONS;

Current_Base::Current_Base(ATOOLS::Flavour* flavs, int n, int* decayindices, string name) :
  Spin_Currents(flavs,n,decayindices), m_n(n), m_name(name), m_flavs(flavs),
  p_i(decayindices)
{
  p_masses = new double[n];
  for(int i=0; i<m_n; ++i) p_masses[i]  = m_flavs[p_i[i]].HadMass();
  msg_Tracking()<<"  Initialized "<<m_name<<" current with "
		<<size()<<" spin combinations"<<endl;
  for(int i=0; i<m_n; i++) {
    msg_Debugging()<<"    flavs["<<i<<"]="<<m_flavs[p_i[i]]<<endl;
    msg_Debugging()<<"    i["<<i<<"]="<<p_i[i]<<endl;
  }
}


Current_Base::~Current_Base()
{
  if (p_i!=NULL) delete[] p_i; p_i=NULL;
  if (p_masses!=NULL) delete[] p_masses; p_masses=NULL;
}

namespace HADRONS {
  void ContractCurrents(HADRONS::Current_Base* c1,
                        HADRONS::Current_Base* c2,
                        const ATOOLS::Vec4D * moms,
                        Complex mefactor,
                        METOOLS::Spin_Amplitudes* amps)
  {
    c1->Calc(moms);
    c2->Calc(moms);
    
    std::vector<int> spins,spins1,spins2;
    for(size_t i=0;i<amps->size();i++) {
      spins=amps->GetSpinCombination(i);
      spins1.clear(); spins2.clear();
      for(size_t j=0;j<c1->GetN();j++)
        spins1.push_back(spins[c1->DecayIndices()[j]]);
      for(size_t j=0;j<c2->GetN();j++)
        spins2.push_back(spins[c2->DecayIndices()[j]]);
      // now we know the spin combinations in both currents
      // let's fill the results:
      (*amps)[i]+=mefactor*c1->Get(spins1)*c2->Get(spins2);
    }
  }

  std::ostream& operator<<(std::ostream& s, const HADRONS::Current_Base& cb)
  {
    s<<cb.Name()<<" current with "<<cb.size()<<" spin combinations:"<<endl;
    for(int i=0; i<cb.size(); i++) {
      s<<"  "<<cb[i]<<endl;
    }
    return s;
  }
}
