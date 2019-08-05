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

/*bool Kernels::SelectOne(Splitting & split,const Mass_Selector * ms) {
  double disc = m_sum * ran->Get(), term;
  for (size_t i=0;i<size();i++) {
    disc -= term = (*this)[i]->Integral(split,ms);
    msg_Out()<<METHOD<<" disc = "<<disc<<" from sum = "<<m_sum<<", "
	     <<"term = "<<term<<" for "
	     <<(*this)[i]->GetFlavs()[0]<<" -> "
	     <<(*this)[i]->GetFlavs()[1]<<" + "
	     <<(*this)[i]->GetFlavs()[2]<<".\n";
    if (disc<0.) {
      p_selected = (*this)[i];
      return true;
    }
  }
  msg_Out()<<METHOD<<" fails, will exit the run for debugging purposes\n";
  exit(1);
  p_selected = NULL;
  return false;
}
*/

ostream & CFPSHOWER::operator<<(ostream &s,Kernels & kernels) {
  for (size_t i=0;i<kernels.size();i++) {
    s<<"     "<<(*kernels[i]);
  }
  return s;
}
