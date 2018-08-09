#include "CFPSHOWER++/Shower/Kernels.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

Kernels::~Kernels() {
  while (!empty()) {
    delete back();
    pop_back();
  }
}

double Kernels::CalcIntegrals(Splitting & split,const Mass_Selector * ms) {
  m_sum = 0.;
  m_integrals.clear();
  double integral;
  for (size_t i=0;i<size();i++) {
    split.SetKernel((*this)[i]);
    m_sum += integral = (*this)[i]->Integral(split,ms);
    m_integrals.push_back(integral);
  }
  //if (m_sum<1.e-6 || split.T0()<1. || split.Q2()<4.) 
  //msg_Out()<<" *** add ["<<size()<<" kernels] = "<<m_sum
  //	     <<" for spectator "<<split.GetSpectator()->Id()
  //	     <<" and Q2 = "<<split.Q2()<<" -> t0 = "<<split.T0()<<"\n";
  return m_sum;
}

bool Kernels::SelectOne() {
  double disc = m_sum * ran->Get();
  for (size_t i=0;i<size();i++) {
    disc -= m_integrals[i];
    if (disc<0.) {
      p_selected = (*this)[i];
      return true;
    }
  }
  msg_Out()<<METHOD<<" fails.\n";
  exit(1);
  p_selected = NULL;
  return false;
}


ostream & CFPSHOWER::operator<<(ostream &s,Kernels & kernels) {
  for (size_t i=0;i<kernels.size();i++) {
    s<<"     "<<(*kernels[i]);
  }
  return s;
}
