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
  m_integrals.assign(size()+1,0.);
  for (size_t i=0;i<size();i++) {
    m_integrals.back() += m_integrals[i] = (*this)[i]->Integral(split,ms);
  }
  return m_integrals.back();
}

ostream & CFPSHOWER::operator<<(ostream &s,Kernels & kernels) {
  for (size_t i=0;i<kernels.size();i++) {
    s<<"     "<<(*kernels[i]);
  }
  return s;
}
