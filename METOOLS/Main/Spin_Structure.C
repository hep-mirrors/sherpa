#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/Blob.H"
#include <iomanip>
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

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
                                 const std::vector<int>& indices) :
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


bool Spin_Amplitudes::IsNLODecay(){
  return false;
}


double Spin_Amplitudes::get_NLO_ME2(){
  return 0.0;
}


std::string METOOLS::Spin_Amplitudes::getType(){
  // Default type for base; derived classes override (e.g., "R", "V", "I").
  return "LO";
}


void Spin_Amplitudes::setBornAmplitude(Spin_Amplitudes* born) {
  // Default: not doing anything
}


void Spin_Amplitudes::setBornAmplitude(std::map<std::string, std::complex<double>> born) {
  // Default: not doing anything
}


std::map<std::string, std::complex<double>> Spin_Amplitudes::getBornAmplitude() {
  return {};
  // Default: not doing anything
}


const PHASIC::Color_Integrator* Spin_Amplitudes::GetColors() const {
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  return nullptr;
}

void Spin_Amplitudes::SetColors(const vector<int> &ci, const vector<int> &cj) {
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}

void Spin_Amplitudes::FixColor(bool yes) {
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}

double Spin_Amplitudes::GetColourWeight(){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}


