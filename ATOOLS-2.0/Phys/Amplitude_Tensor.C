#include "Amplitude_Tensor.H"
#include "Message.H"
#include <iomanip>
#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace std;

Amplitude_Tensor::Amplitude_Tensor( Particle_Vector particles) :
  m_particles(particles), p_colormatrix(NULL)
{
  size_t n=1;
  for(size_t i=0;i<m_particles.size();i++) {
    Flavour flav = m_particles[i]->Flav();
    int spincombinations = flav.IntSpin()+1;
    n*=spincombinations;
  }
  m_amplitudes = std::vector<std::vector<Complex> >(n);
}

Amplitude_Tensor::~Amplitude_Tensor()
{
}

size_t Amplitude_Tensor::GetAmplitudeNumber(const vector<int>& spins) const
{
  if(spins.size()!=m_particles.size()) {
    msg.Error()<<METHOD<<" Error: wrong size of spin vector."<<endl;
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
    msg.Error()<<METHOD<<" Error: wrong size of spin vector."<<endl;
  }
  int mult(1);
  size_t num(0);
  for(size_t i=spins.size()-1; i+1>0; i--) {
    num += mult * spins[i].second;
    mult *= (m_particles[i]->Flav().IntSpin()+1);
  }
  if(num>Size()) {
    msg.Error()<<METHOD<<" Error: tried to access amplitude out of bounce. "
      <<"num="<<num<<" > "<<Size()<<endl;
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
  
  for(size_t i=0;i<newamps.Size();i++) {
    std::vector<int> spins = newamps.GetSpinCombination(i);
    vector<Complex> amp(amps1->ColorSize()*amps2->ColorSize(),Complex(0.0,0.0));
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

//       amp+=amps1->GetAmplitude(spins1)*amps2->GetAmplitude(spins2);

      vector<Complex> amps1_colors = amps1->GetAmplitude(spins1);
      vector<Complex> amps2_colors = amps2->GetAmplitude(spins2);
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
        msg.Error()<<METHOD<<" Error: contraction not implemented for two "
          <<"colorful Amplitude_Tensors yet."<<endl;
        abort();
      }
    }
    newamps.InsertAmplitude(amp,i);
  }
  return newamps;
}

void Amplitude_Tensor::InsertAmplitude(Complex amp, vector<pair<int,int> >& spins)
{
  /*! only one color combination in that one spin combination*/
  vector<Complex> amps(1,amp);
  m_amplitudes[GetAmplitudeNumber(spins)]=amps;
}

void Amplitude_Tensor::InsertAmplitude(vector<Complex> amps, vector<pair<int,int> >& spins)
{
  /*! amplitudes for multiple color combinations in that one spin combination*/
  m_amplitudes[GetAmplitudeNumber(spins)]=amps;
}

void Amplitude_Tensor::InsertAmplitude(Complex amp, size_t index)
{
  /*! only one color combination in that one spin combination */
  vector<Complex> amps(1,amp);
  m_amplitudes[index]=amps;
}

void Amplitude_Tensor::InsertAmplitude(vector<Complex> amps, size_t index)
{
  /*! amplitudes for multiple color combinations in that one spin combination*/
  m_amplitudes[index]=amps;
}

vector<Complex> Amplitude_Tensor::GetAmplitude(const std::vector<int>& spins) const
{
  return m_amplitudes[GetAmplitudeNumber(spins)];
}

vector<Complex> Amplitude_Tensor::GetAmplitude(size_t index) const
{
  return m_amplitudes[index];
}

double Amplitude_Tensor::SumSquare() const
{
  Complex value(0.0,0.0);
  for( size_t i=0; i<Size(); ++i ) {
    for( size_t j=0; j<ColorSize(); j++) {
      for( size_t k=0; k<ColorSize(); k++) {
        Complex Mj = m_amplitudes[i][j];
        Complex Mk = conj(m_amplitudes[i][k]);
        Complex Mjk = (*p_colormatrix)[j][k];
        value += Mj*Mk*Mjk; // should be real?
//         if(!IsZero((value.imag()))) PRINT_INFO("value.imag()="<<value.imag());
      }
    }
  }
  return value.real();
}

void Amplitude_Tensor::Recreate( Amplitude_Tensor* newamps)
{
  m_particles = newamps->m_particles;
  m_amplitudes = newamps->m_amplitudes;
  p_colormatrix = newamps->p_colormatrix;
}

void Amplitude_Tensor::CreateTrivial()
{
  size_t n = Size();
  if(n<1) {
    msg.Error()<<METHOD<<" Error: m_amplitudes has size "<<n<<endl;
    abort();
  }
  m_amplitudes = vector<vector<Complex> >(n,vector<Complex>(1,Complex(1.0,0.0)));
}

size_t Amplitude_Tensor::Size() const
{
  return m_amplitudes.size();
}

size_t Amplitude_Tensor::ColorSize() const
{
  return p_colormatrix->Rank();
}

CMatrix* Amplitude_Tensor::GetColorMatrix() const
{
  return p_colormatrix;
}

void Amplitude_Tensor::SetColorMatrix(CMatrix* colormatrix)
{
  p_colormatrix = colormatrix;
}

std::ostream& ATOOLS::operator<<( std::ostream& ostr, const Amplitude_Tensor & amps) {
  ostr<<"Amplitude_Tensor with "<<amps.m_particles.size()<<" particles and "
    <<amps.m_amplitudes.size()<<" spin combinations:"<<endl;
  for(size_t i=0;i<amps.m_particles.size();i++) {
    ostr<<setw(8)<<amps.m_particles[i]->Flav()<<" | ";
  }
  ostr<<"first color amplitude"<<endl;
  for(size_t i=0;i<amps.m_amplitudes.size();i++) {
    std::vector<int> spins = amps.GetSpinCombination(i);
    for(size_t j=0;j<amps.m_particles.size();j++) {
      ostr<<setw(8)<<spins[j]<<" | ";
    }
    ostr<<amps.m_amplitudes[i][0]<<endl;
  }
  return ostr;
}

namespace ATOOLS {
  template <> Blob_Data<Amplitude_Tensor*>::~Blob_Data() {
    if(m_data) delete m_data; m_data=NULL;
  }
  template class Blob_Data<Amplitude_Tensor*>;
  template Amplitude_Tensor* &Blob_Data_Base::Get<Amplitude_Tensor*>();
}
