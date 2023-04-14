#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"

//#define COMPILE__Getter_Function
//#define OBJECT_TYPE NEUTRINOS::Scatter_Current_Base
//#define PARAMETER_TYPE NEUTRINOS::ME_Parameters
#include "ATOOLS/Org/Getter_Function.C"

using namespace NEUTRINOS;
using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Scatter_Current_Base::Scatter_Current_Base(const ATOOLS::Flavour_Vector& flavs,
					   const std::vector<int>& indices,
					   const std::string& name) :
  Spin_Structure<ATOOLS::Vec4C>(flavs,indices),
  m_flavs(flavs), 
  m_name(name)
{ //This is the old unedited code.
  m_indices.resize(indices.size());
  m_masses.resize(indices.size());
  for(int i=0; i<indices.size(); ++i) {
    m_indices[i] = indices[i];
    m_masses[i]  = m_flavs[m_indices[i]].HadMass();
  }
  msg_Info()<<METHOD<<" initialized "<<m_name<<" current with "<<size()<<" spin combinations:\n";
  for(int i=0; i<m_indices.size(); i++) 
    msg_Info()<<"    flavs["<<i<<" --> "<<m_indices[i]<<"] = "<<m_flavs[m_indices[i]]<<"\n";
}
