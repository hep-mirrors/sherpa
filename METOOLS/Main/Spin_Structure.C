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


const PHASIC::Color_Integrator* Spin_Amplitudes::GetColorIntegrator() const {
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
  return 0.0;
}

ATOOLS::Vec4D_Vector Spin_Amplitudes::GetMomenta() {
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  return ATOOLS::Vec4D_Vector(0);
}


double Spin_Amplitudes::getColourFactor(const PHASIC::Color_Integrator* p_ci){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  return 1.0;
}


double Spin_Amplitudes::getSign(){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  return 1.0;
}


std::array<int, 3> Spin_Amplitudes::getDipoleIndices(){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  std::array<int, 3> a = {0, 0, 0};
  return a;
}


void Spin_Amplitudes::MergeDiagrams(const METOOLS::Spin_Amplitudes* second_diagram){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}


void Spin_Amplitudes::SetFullME2(double fullME2){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
}


double Spin_Amplitudes::GetFullME2(){
  msg_Error()<<METHOD<<": Virtual function called."<<endl;
  Abort();
  return 0.0;
}

