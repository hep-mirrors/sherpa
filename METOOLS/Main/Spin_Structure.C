#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/Blob.H"
#include <iomanip>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

namespace METOOLS {
bool SortByFirst(const pair<int,int> p1, const pair<int,int> p2) {
  return p1.first < p2.first;
}
}

Spin_Amplitudes::~Spin_Amplitudes() {}

Spin_Amplitudes::Spin_Amplitudes(const std::vector<int>& spins,
                                 const Complex& value) :
  Spin_Structure<Complex>(spins, value)
{
}

Spin_Amplitudes::Spin_Amplitudes(const Particle_Vector& particles) :
  Spin_Structure<Complex>(particles)
{
}

Spin_Amplitudes::Spin_Amplitudes(const Flavour_Vector& flavs,
                                 const Complex& value) :
  Spin_Structure<Complex>(flavs,value)
{
}

Spin_Amplitudes::Spin_Amplitudes(const Flavour_Vector& flavs,
                                 const vector<int>& indices) :
  Spin_Structure<Complex>(flavs,indices)
{
}

double Spin_Amplitudes::SumSquare() const {
  double value(0);
  for( size_t i=0; i<this->size(); i++ ) {
    value += norm( (*this)[i] );
  }
  return value;
}

void Spin_Amplitudes::Calculate(const Vec4D_Vector& momenta, bool anti) {
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}
