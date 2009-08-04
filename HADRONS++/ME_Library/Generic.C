#include "HADRONS++/ME_Library/Generic.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

#include "METOOLS/Main/Partial_Amplitude_Base.H"

Generic::Generic(ATOOLS::Flavour* flavs, int n, int* decayindices, 
                 std::string name):
  HD_ME_Base(flavs,n,decayindices,name)
{
  p_me=Partial_Amplitude_Base::Select(flavs,n);
}

Generic::~Generic() {
  delete p_me;
}

void Generic::operator()(const ATOOLS::Vec4D* p, METOOLS::Spin_Amplitudes* amps)
{
  (*p_me)(p, m_anti);
  if(amps->size()!=p_me->size())
    THROW(fatal_error, "amps->size()!=p_me->size()");
  for(size_t i(0);i<amps->size();++i) {
    (*amps)[i]=(*p_me)[i];
  }
}

DEFINE_ME_GETTER(Generic,Generic_Getter,"Generic")

void Generic_Getter::PrintInfo(std::ostream &str,const size_t width) const {
  str<<"Chooses a generic matrix element according to the spin structure.";
}
