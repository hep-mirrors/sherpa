#include "Amplitude_Tensor.H"
#include "Message.H"
#include <iomanip>
#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace std;

Amplitude_Tensor::Amplitude_Tensor( Particle_Vector particles) :
  m_particles(particles)
{
  size_t n=1;
  for(size_t i=0;i<m_particles.size();i++) {
    Flavour flav = m_particles[i]->Flav();
    int spincombinations = flav.IntSpin()+1;
    n*=spincombinations;
  }
  m_amplitudes = std::vector<Complex>(n);
}

size_t Amplitude_Tensor::GetAmplitudeNumber(const vector<int>& spins) const
{
  if(spins.size()!=m_particles.size()) {
    msg.Error()<<METHOD<<" Error: wrong size of spin vector:"<<endl;
    msg.Error()<<"spins.size()="<<spins.size()<<endl;
    for(size_t i=0;i<spins.size();i++) {
      msg.Error()<<"spins["<<i<<"]="<<spins[i]<<endl;
    }
    msg.Error()<<"m_particles.size()="<<m_particles.size()<<endl;
    for(size_t i=0;i<m_particles.size();i++) {
      msg.Error()<<"m_particles["<<i<<"]="<<m_particles[i]->Flav()<<endl;
    }
    abort();
  }
  int mult(1);
  size_t num(0);
  for(size_t i=spins.size()-1; i+1>0; i--) {
    num += mult * spins[i];
    mult *= (m_particles[i]->Flav().IntSpin()+1);
  }
  return num;
}

bool SortByFirst(const pair<int,int> p1, const pair<int,int> p2) {
  return p1.first < p2.first;
}

size_t Amplitude_Tensor::GetAmplitudeNumber(vector<pair<int,int> >& spins) const
{
  sort(spins.begin(),spins.end(),SortByFirst);
  
  if(spins.size()!=m_particles.size()) {
    msg.Error()<<METHOD<<" Error: wrong size of spin vector:"<<endl;
    msg.Error()<<"spins.size()="<<spins.size()<<endl;
    for(size_t i=0;i<spins.size();i++) {
      msg.Error()<<"spins["<<i<<"]="<<spins[i].second<<endl;
    }
    msg.Error()<<"m_particles.size()="<<m_particles.size()<<endl;
    for(size_t i=0;i<m_particles.size();i++) {
      msg.Error()<<"m_particles["<<i<<"]="<<m_particles[i]->Flav()<<endl;
    }
    abort();
  }
  int mult(1);
  size_t num(0);
  for(size_t i=spins.size()-1; i+1>0; i--) {
    num += mult * spins[i].second;
    mult *= (m_particles[i]->Flav().IntSpin()+1);
  }
  if(num>m_amplitudes.size()) {
    msg.Error()<<METHOD<<" Error: tried to access amplitude out of bounce. "
      <<"num="<<num<<" > "<<m_amplitudes.size()<<endl;
  }
  return num;
}

std::vector<int> Amplitude_Tensor::GetSpinCombination(size_t ampnumber) const
{
  std::vector<int> spins(m_particles.size());
  for(int i=m_particles.size()-1;i>=0;--i) {
    int spincombinations = m_particles[i]->Flav().IntSpin()+1;
    spins[i] = ampnumber%spincombinations;
    ampnumber = (ampnumber-spins[i])/spincombinations;
  }
  return spins;
}

Amplitude_Tensor ATOOLS::Contraction(Particle* part,
                                     Amplitude_Tensor* const amps1,
                                     Amplitude_Tensor* const amps2
                                    )
{
  Particle_Vector remainingparts;
  for(size_t i=0;i<amps1->m_particles.size();i++) {
    if(part!=amps1->m_particles[i]) remainingparts.push_back(amps1->m_particles[i]);
  }
  for(size_t i=0;i<amps2->m_particles.size();i++) {
    if(part!=amps2->m_particles[i]) remainingparts.push_back(amps2->m_particles[i]);
  }
  if(remainingparts.size()!=amps1->m_particles.size()-1+amps2->m_particles.size()-1) {
    msg.Error()<<METHOD<<" Error: not exactly one particle to contract in event "
      <<rpa.gen.NumberOfDicedEvents()<<"."<<endl;
    for(size_t i=0;i<remainingparts.size();i++) {
      msg.Error()<<"remainingparts["<<i<<"]="<<remainingparts[i]->Flav()<<endl;
    }
    for(size_t i=0;i<amps1->m_particles.size();i++) {
      msg.Error()<<"amps1->m_particles["<<i<<"]="<<amps1->m_particles[i]->Flav()<<endl;
    }
    for(size_t i=0;i<amps2->m_particles.size();i++) {
      msg.Error()<<"amps2->m_particles["<<i<<"]="<<amps2->m_particles[i]->Flav()<<endl;
    }
    msg.Error()<<"contracting particle:"<<endl<<(*part)<<endl;
//     abort();
  }
  Amplitude_Tensor newamps(remainingparts);
  
  for(size_t i=0;i<newamps.m_amplitudes.size();i++) {
    std::vector<int> spins = newamps.GetSpinCombination(i);
    Complex amp(0.0,0.0);
    for(int j=0;j<part->Flav().IntSpin()+1;j++) {
      
      vector<int> spins1;
      int offset=0;
      for(size_t k=0;k<amps1->m_particles.size();k++) {
        if(amps1->m_particles[k]==part) {
          spins1.push_back(j);
          offset=1;
        }
        else spins1.push_back(spins[k-offset]);
      }
      vector<int> spins2;
      offset=0;
      for(size_t k=0;k<amps2->m_particles.size();k++) {
        if(amps2->m_particles[k]==part) {
          spins2.push_back(j);
          offset=1;
        }
        else spins2.push_back(spins[amps1->m_particles.size()-1+k-offset]);
      }
      amp+=amps1->GetAmplitude(spins1)*amps2->GetAmplitude(spins2);
    }
    newamps.InsertAmplitude(amp,i);
  }
  return newamps;
}

void Amplitude_Tensor::InsertAmplitude(Complex amp, const vector<int>& spins)
{
  m_amplitudes[GetAmplitudeNumber(spins)]=amp;
}

void Amplitude_Tensor::InsertAmplitude(Complex amp, vector<pair<int,int> >& spins)
{
  m_amplitudes[GetAmplitudeNumber(spins)]=amp;
}

void Amplitude_Tensor::InsertAmplitude(Complex amp, size_t index)
{
  m_amplitudes[index]=amp;
}

Complex Amplitude_Tensor::GetAmplitude(const std::vector<int>& spins) const
{
  return m_amplitudes[GetAmplitudeNumber(spins)];
}

Complex Amplitude_Tensor::GetAmplitude(size_t index) const
{
  return m_amplitudes[index];
}

int Amplitude_Tensor::Size()
{
  return m_amplitudes.size();
}

double Amplitude_Tensor::SumSquare() const
{
  double value=0.0;
  for( size_t i=0; i<m_amplitudes.size(); ++i ) {
    value += norm(m_amplitudes[i]);
  }
  return value;
}

void Amplitude_Tensor::Recreate( Amplitude_Tensor* newamps)
{
  m_particles = newamps->m_particles;
  m_amplitudes = newamps->m_amplitudes;
}

void Amplitude_Tensor::CreateTrivial()
{
  size_t n = m_amplitudes.size();
  if(n<1) {
    msg.Error()<<METHOD<<" Error: m_amplitudes has size "<<n<<endl;
    abort();
  }
  m_amplitudes = std::vector<Complex>(n,Complex(1.0,0.0));
}

std::ostream& ATOOLS::operator<<( std::ostream& ostr, const Amplitude_Tensor & amps) {
  ostr<<"Amplitude_Tensor with "<<amps.m_particles.size()<<" particles and "
    <<amps.m_amplitudes.size()<<" spin combinations:"<<endl;
  for(size_t i=0;i<amps.m_particles.size();i++) {
    ostr<<setw(8)<<amps.m_particles[i]->Flav()<<" | ";
  }
  ostr<<"amplitudes"<<endl;
  for(size_t i=0;i<amps.m_amplitudes.size();i++) {
    std::vector<int> spins = amps.GetSpinCombination(i);
    for(size_t j=0;j<amps.m_particles.size();j++) {
      ostr<<setw(8)<<spins[j]<<" | ";
    }
    ostr<<amps.m_amplitudes[i]<<endl;
  }
  return ostr;
}

namespace ATOOLS {
  template <> Blob_Data<Amplitude_Tensor*>::~Blob_Data() { delete m_data; }
  template class Blob_Data<Amplitude_Tensor*>;
  template Amplitude_Tensor* &Blob_Data_Base::Get<Amplitude_Tensor*>();
}
