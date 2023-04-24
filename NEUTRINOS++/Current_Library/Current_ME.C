#include "NEUTRINOS++/Current_Library/Current_ME.H"
#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "METOOLS/Main/XYZFuncs.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

Current_ME::Current_ME(const ATOOLS::Flavour_Vector& flavs,
                       const std::vector<int>& decayindices,
                       const std::string& name):
  HADRONS::HD_ME_Base(flavs,decayindices,name), p_c1(NULL), p_c2(NULL), m_factor(1.) {
  ffs->Output();
}

Current_ME::~Current_ME() {
  if (p_c1) delete p_c1;
  if (p_c2) delete p_c2;
}

void Current_ME::Calculate(const Vec4D_Vector& momenta,METOOLS::XYZFunc * F)
{
  p_c1->Calc(momenta, F);
  p_c2->Calc(momenta, F);

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

void Current_ME::Calculate(const Vec4D_Vector& momenta,bool m_anti)
{
  //JW: Duplicate momenta
  Vec4D_Vector momenta_doubled = momenta;
  momenta_doubled.insert(momenta_doubled.end(), momenta.begin(), momenta.end());
  const Vec4D_Vector momenta_xyz = momenta_doubled;

  vector<int> p_i_doubled = p_i;
  p_i_doubled.insert(p_i_doubled.end(), p_i.begin(), p_i.end());
  for(int i=p_i.size(); i<2*p_i.size(); i++) {
    //Index new anti-particles with a shift +n for n original particles
    p_i_doubled[i] += p_i.size();
  }
  const vector<int> p_i_xyz = p_i_doubled;

  //JW: Duplicate flavours
  ATOOLS::Flavour_Vector flavs_doubled = m_flavs;
  flavs_doubled.insert(flavs_doubled.end(), m_flavs.begin(), m_flavs.end());
  for(int i=m_flavs.size(); i<2*m_flavs.size(); i++) {
    //Replace second half of flavours with anti-particle version.
    flavs_doubled[i] = flavs_doubled[i].Bar();
  }
  const ATOOLS::Flavour_Vector m_flavs_xyz = flavs_doubled;

  METOOLS::XYZFunc F(momenta_xyz, m_flavs_xyz, false, p_i_xyz);
  Calculate(momenta_xyz, &F);
}

double Current_ME::operator()(const Vec4D_Vector& momenta) {
  Calculate(momenta);
  return SumSquare();
}
