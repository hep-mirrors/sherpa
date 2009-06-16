#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::HD_ME_Base
#define PARAMETER_TYPE HADRONS::Flavour_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

HD_ME_Base::HD_ME_Base(Flavour * flavs, int n, int* decayindices, string name) :
  m_n(n), m_name(name), m_flavs(flavs), p_i(decayindices)
{
  p_masses = new double[n];
  p_masses2 = new double[n];
  
  for(int meindex=0; meindex<m_n; meindex++) {
    p_masses[meindex]  = m_flavs[p_i[meindex]].HadMass();
    p_masses2[meindex] = p_masses[meindex]*p_masses[meindex];
  }
  msg_Tracking()<<"  Initialized "<<m_name<<" ME."<<endl;
  for(int i=0; i<m_n; i++) {
    msg_Debugging()<<"    flavs["<<i<<"]="<<m_flavs[i]<<endl;
    msg_Debugging()<<"    i["<<i<<"]="<<p_i[i]<<endl;
  }
}

HD_ME_Base::~HD_ME_Base()
{
  if (p_masses)  { delete [] p_masses;  p_masses  = NULL; }
  if (p_masses2) { delete [] p_masses2; p_masses2 = NULL; }
  if (p_i!=NULL) { delete[] p_i; p_i=NULL; }
}

double HD_ME_Base::lambdaNorm(const double M,const double m1,const double m2) {
  return sqrt((sqr(M)-sqr(m1+m2))*(sqr(M)-sqr(m1-m2)))/(2.*M);
}

bool HD_ME_Base::SetColorFlow(std::vector<Particle*> outparts,int n_q, int n_g)
{
  return false;
}
