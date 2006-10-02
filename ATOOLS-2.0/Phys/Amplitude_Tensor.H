#ifndef Amplitude_Tensor_h
#define Amplitude_Tensor_h

#include "MyComplex.H"
#include "Particle.H"
#include "Blob.H"

namespace ATOOLS
{
  class Amplitude_Tensor {
    std::vector<Complex>  m_amplitudes;
    Particle_Vector       m_particles;
    
    size_t GetAmplitudeNumber(const std::vector<int>& spins) const;
    size_t GetAmplitudeNumber(std::vector<std::pair<int,int> >& spins) const;
    std::vector<int> GetSpinCombination(size_t ampnumber) const;
  public:
    Amplitude_Tensor(std::vector<Particle* > particles);
    
    friend Amplitude_Tensor Contraction(
                                        Particle* part,
                                        Amplitude_Tensor* const amps1,
                                        Amplitude_Tensor* const amps2
                                       );

    void    InsertAmplitude(Complex amp, const std::vector<int>& spins);
    void    InsertAmplitude(Complex amp, std::vector<std::pair<int,int> >& spins);
    void    InsertAmplitude(Complex amp, size_t index);
    Complex GetAmplitude(const std::vector<int>& spins) const;
    Complex GetAmplitude(size_t index) const;
    double  SumSquare() const;
    void    Recreate(Amplitude_Tensor* newamps);
    void    CreateTrivial();
    int     Size();
    friend std::ostream& operator<<( std::ostream&, const Amplitude_Tensor &);
  };
}
#endif
