#include "METOOLS/HadronCurrents/F0_Decays.H"
#include "METOOLS/HadronCurrents/Scalar_Decays.H"
#include "ATOOLS/Phys/Flavour.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

F_0_500_Lineshape::F_0_500_Lineshape() :
  Total_Width_Base(Flavour(kf_f_0_600)) {
  vector<Flavour> outflavs;
  outflavs = { Flavour(kf_pi_plus), Flavour(kf_pi_plus).Bar() };
  Partial_Width_Base * f02pipiC = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(f02pipiC);
  outflavs = { Flavour(kf_pi), Flavour(kf_pi) };
  Partial_Width_Base * f02pipiN = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(f02pipiN);
}
