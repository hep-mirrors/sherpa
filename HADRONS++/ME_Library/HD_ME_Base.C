#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::HD_ME_Base
#define PARAMETER_TYPE HADRONS::ME_Parameters
#include "ATOOLS/Org/Getter_Function.C"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

HD_ME_Base::HD_ME_Base(const ATOOLS::Flavour_Vector& flavs,
                       const std::vector<int>& decayindices,
                       const std::string& name) :
  METOOLS::Spin_Amplitudes(flavs,Complex(0.0,0.0)),
  m_name(name), m_flavs(flavs), m_factor(1.0)
{
  p_masses = new double[m_flavs.size()];
  p_masses2 = new double[m_flavs.size()];

  p_i.resize(m_flavs.size());
  for(int meindex=0; meindex<m_flavs.size(); meindex++) {
    p_i[meindex]=decayindices[meindex];
    p_masses[meindex]  = m_flavs[p_i[meindex]].HadMass();
    p_masses2[meindex] = p_masses[meindex]*p_masses[meindex];
  }
  msg_Tracking()<<"  Initialized "<<m_name<<" ME."<<endl;
  for(int i=0; i<m_flavs.size(); i++) {
    msg_Debugging()<<"    flavs["<<i<<"]="<<m_flavs[i]<<endl;
    msg_Debugging()<<"    i["<<i<<"]="<<p_i[i]<<endl;
  }
}

HD_ME_Base::~HD_ME_Base()
{
  if (p_masses)  { delete [] p_masses;  p_masses  = NULL; }
  if (p_masses2) { delete [] p_masses2; p_masses2 = NULL; }
}

double HD_ME_Base::lambdaNorm(const double M,const double m1,const double m2) {
  return sqrt((sqr(M)-sqr(m1+m2))*(sqr(M)-sqr(m1-m2)))/(2.*M);
}

bool HD_ME_Base::SetColorFlow(std::vector<Particle*> outparts,int n_q, int n_g,
                              bool anti)
{
  return false;
}
