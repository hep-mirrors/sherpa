#include "HELICITIES/Main/Spin_Structure.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include <iomanip>

using namespace HELICITIES;
using namespace ATOOLS;
using namespace std;

namespace HELICITIES {
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

Spin_Amplitudes::Spin_Amplitudes(const Flavour* flavs, size_t size,
                                 const Complex& value) :
  Spin_Structure<Complex>(flavs,size,value)
{
}

Spin_Amplitudes::Spin_Amplitudes(const Flavour_Vector& flavs,
                                 const Complex& value) :
  Spin_Structure<Complex>(flavs,value)
{
}

double Spin_Amplitudes::SumSquare() const {
  double value(0);
  for( size_t i=0; i<this->size(); i++ ) {
    value += norm( (*this)[i] );
  }
  return value;
}

void Spin_Amplitudes::Calculate(const ATOOLS::Vec4D* momenta, bool anti) {
  THROW(fatal_error, "Virtual function called.");
}



Amplitude_Tensor::Amplitude_Tensor( Particle_Vector particles) :
  Spin_Structure<std::vector<Complex> >(particles), m_particles(particles),
  p_colormatrix(NULL), p_contracted(NULL)
{
}

Amplitude_Tensor::~Amplitude_Tensor()
{
  if(p_contracted) { delete p_contracted; p_contracted=NULL; }
}

namespace HELICITIES {
Amplitude_Tensor Contraction(Particle* part1, Particle* part2,
                                     Amplitude_Tensor* const amps1,
                                     Amplitude_Tensor* const amps2
                                    )
{
  if(part1->Flav()!=part2->Flav()) {
    msg_Error()<<METHOD<<" Contracting particles are not the same:"<<endl
      <<(*part1)<<endl
      <<(*part2)<<endl;
  }
  Particle_Vector remainingparts;
  for(size_t i=0;i<amps1->m_particles.size();i++) {
    if(part1!=amps1->m_particles[i]) remainingparts.push_back(amps1->m_particles[i]);
  }
  for(size_t i=0;i<amps2->m_particles.size();i++) {
    if(part2!=amps2->m_particles[i]) remainingparts.push_back(amps2->m_particles[i]);
  }
  if(remainingparts.size()!=amps1->m_particles.size()-1+amps2->m_particles.size()-1) {
    msg_Error()<<METHOD<<" Error: not exactly one particle to contract."<<endl;
    for(size_t i=0;i<remainingparts.size();i++) {
      msg_Error()<<"remainingparts["<<i<<"]="<<remainingparts[i]->Flav()<<endl;
    }
    for(size_t i=0;i<amps1->m_particles.size();i++) {
      msg_Error()<<"amps1->m_particles["<<i<<"]="<<amps1->m_particles[i]->Flav()<<endl;
    }
    for(size_t i=0;i<amps2->m_particles.size();i++) {
      msg_Error()<<"amps2->m_particles["<<i<<"]="<<amps2->m_particles[i]->Flav()<<endl;
    }
    msg_Error()<<"contracting particle1:"<<endl<<(*part1)<<endl
      <<"contracting particle2:"<<endl<<(*part2)<<endl;
//     abort();
  }
  Amplitude_Tensor newamps(remainingparts);

  for(size_t i=0;i<newamps.size();i++) {
    std::vector<int> spins = newamps.GetSpinCombination(i);
    vector<Complex> amp(amps1->ColorSize()*amps2->ColorSize(),Complex(0.0,0.0));
    for(int j=0;j<part1->Flav().IntSpin()+1;j++) {

      vector<int> spins1;
      int offset=0;
      for(size_t k=0;k<amps1->m_particles.size();k++) {
        if(amps1->m_particles[k]==part1) {
          spins1.push_back(j);
          offset=1;
        }
        else spins1.push_back(spins[k-offset]);
      }
      vector<int> spins2;
      offset=0;
      for(size_t k=0;k<amps2->m_particles.size();k++) {
        if(amps2->m_particles[k]==part2) {
          spins2.push_back(j);
          offset=1;
        }
        else spins2.push_back(spins[amps1->m_particles.size()-1+k-offset]);
      }

//       amp+=amps1->GetAmplitude(spins1)*amps2->GetAmplitude(spins2);

      vector<Complex> amps1_colors = amps1->Get(spins1);
      vector<Complex> amps2_colors = amps2->Get(spins2);
      size_t m=0;
      for(size_t k=0;k<amps1_colors.size();k++) {
        for(size_t l=0;l<amps2_colors.size();l++) {
          amp[m] += amps1_colors[k]*amps2_colors[l];
          m++;
        }
      }
      // todo: more general for both ColorSize() > 1
      if(amps1->ColorSize()==1) {
        newamps.SetColorMatrix(amps2->GetColorMatrix());
      }
      else if(amps2->ColorSize()==1) {
        newamps.SetColorMatrix(amps1->GetColorMatrix());
      }
      else {
        msg_Error()<<METHOD<<" Error: contraction not implemented for two "
          <<"colorful Amplitude_Tensors yet."<<endl;
        abort();
      }
    }
    newamps.Insert(amp,i);
  }
  return newamps;
}
}

void Amplitude_Tensor::Contract( Particle* part1, Amplitude_Tensor* const amps, Particle* part2 )
{
  if(p_contracted) {
    Amplitude_Tensor* todelete = p_contracted;
    (*this)=(*p_contracted);
    delete todelete;
    p_contracted=NULL;
  }
  else (*this) = Contraction(part1,part2,this,amps);
}

double Amplitude_Tensor::SoftContract( Particle* part1,
                                       Amplitude_Tensor* const amps,
                                       Particle* part2 )
{
  if(p_contracted)
    (*p_contracted) = Contraction(part1,part2,this,amps);
  else
    p_contracted = new Amplitude_Tensor(Contraction(part1,part2,this,amps));
  return p_contracted->SumSquare();
}

double Amplitude_Tensor::SumSquare() const
{
  Complex value(0.0,0.0);
  for( size_t i=0; i<this->size(); i++ ) {
    for( size_t j=0; j<ColorSize(); j++) {
      for( size_t k=0; k<ColorSize(); k++) {
        Complex Mj = (*this)[i][j];
        Complex Mk = conj((*this)[i][k]);
        Complex Mjk = (*p_colormatrix)[j][k];
        value += Mj*Mk*Mjk;
      }
    }
  }
  // should be real
  if(!IsZero((value.imag()))) PRINT_INFO("value.imag()="<<value.imag());
  return value.real();
}

size_t Amplitude_Tensor::ColorSize() const
{
  if(!p_colormatrix) {
    msg_Error()<<METHOD<<" Warning: asking for ColorSize(), but p_colormatrix has "
      <<"not been set yet."<<endl;
    return 1;
  }
  return p_colormatrix->Rank();
}

const CMatrix* Amplitude_Tensor::GetColorMatrix() const
{
  return p_colormatrix;
}

void Amplitude_Tensor::SetColorMatrix(const CMatrix* colormatrix)
{
  p_colormatrix = colormatrix;
}

bool Amplitude_Tensor::Contains(const Particle* part) const
{
  Particle_Vector::const_iterator it;
  for(it=m_particles.begin();it!=m_particles.end();it++) {
    if( part == (*it) ) return true;
  }
  return false;
}

const Particle_Vector& Amplitude_Tensor::Particles() const
{
  return m_particles;
}

namespace HELICITIES {
std::ostream& operator<<( std::ostream& ostr, const Amplitude_Tensor & amps) {
  ostr<<"   Amplitude_Tensor with "<<amps.m_particles.size()<<" particles and "
    <<amps.size()<<" spin combinations:"<<endl;
  ostr<<"   ";
  for(size_t i=0;i<amps.m_particles.size();i++) {
    ostr<<setw(8)<<amps.m_particles[i]->Flav()<<" | ";
  }
  ostr<<"first color amplitude"<<endl;
  for(size_t i=0;i<amps.size();i++) {
    ostr<<setw(3)<<i;
    std::vector<int> spins = amps.GetSpinCombination(i);
    for(size_t j=0;j<amps.m_particles.size();j++) {
      ostr<<setw(8)<<spins[j]<<" | ";
    }
    ostr<<amps[i][0]<<endl;
  }
  return ostr;
}
}

namespace ATOOLS {
  template <> Blob_Data<Amplitude_Tensor*>::~Blob_Data() {
    if(m_data) delete m_data; m_data=NULL;
  }
  template class Blob_Data<Amplitude_Tensor*>;
  template Amplitude_Tensor* &Blob_Data_Base::Get<Amplitude_Tensor*>();
}
