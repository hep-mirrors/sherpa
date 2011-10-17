#include "HADRONS++/Current_Library/Current_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::Current_Base
#define PARAMETER_TYPE HADRONS::ME_Parameters
#include "ATOOLS/Org/Getter_Function.C"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;
using namespace HADRONS;

Current_Base::Current_Base(const ATOOLS::Flavour_Vector& flavs,
                           const std::vector<int>& decayindices,
                           const std::string& name) :
  Spin_Structure<ATOOLS::Vec4C>(flavs,decayindices),
  m_flavs(flavs),
  m_name(name)
{
  p_i.resize(decayindices.size());
  p_masses = new double[decayindices.size()];
  for(int i=0; i<decayindices.size(); ++i) {
    p_i[i] = decayindices[i];
    p_masses[i] = m_flavs[p_i[i]].HadMass();
  }
  msg_Tracking()<<"  Initialized "<<m_name<<" current with "
		<<size()<<" spin combinations"<<endl;
  for(int i=0; i<p_i.size(); i++) {
    msg_Debugging()<<"    flavs["<<i<<"]="<<m_flavs[p_i[i]]<<endl;
    msg_Debugging()<<"    i["<<i<<"]="<<p_i[i]<<endl;
  }
}


Current_Base::~Current_Base()
{
  if (p_masses!=NULL) delete[] p_masses; p_masses=NULL;
}

namespace HADRONS {
  std::ostream& operator<<(std::ostream& s, const HADRONS::Current_Base& cb)
  {
    s<<cb.Name()<<" current with "<<cb.size()<<" spin combinations:"<<endl;
    for(size_t i=0; i<cb.size(); i++) {
      s<<"  "<<cb[i]<<endl;
    }
    return s;
  }
}
