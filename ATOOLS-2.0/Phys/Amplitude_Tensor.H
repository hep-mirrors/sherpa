#ifndef Amplitude_Tensor_h
#define Amplitude_Tensor_h

#include "MyComplex.H"
#include "Particle.H"
#include "Blob.H"
#include "Matrix.H"

namespace ATOOLS
{
  class Amplitude_Tensor {
    std::vector<std::vector<Complex> >  m_amplitudes;
    Particle_Vector                     m_particles;
    const CMatrix                     * p_colormatrix;
    
    size_t GetAmplitudeNumber(const std::vector<int>& spins) const;
    size_t GetAmplitudeNumber(std::vector<std::pair<int,int> >& spins) const;
    std::vector<int> GetSpinCombination(size_t ampnumber) const;
  public:
    Amplitude_Tensor(Particle_Vector particles);
    ~Amplitude_Tensor();
    
    friend Amplitude_Tensor Contraction(
                                        Particle* part1, Particle* part2,
                                        Amplitude_Tensor* const amps1,
                                        Amplitude_Tensor* const amps2
                                       );

    void    InsertAmplitude(Complex amp, std::vector<std::pair<int,int> >& spins);
    void    InsertAmplitude(std::vector<Complex> amp, std::vector<std::pair<int,int> >& spins);
    
    void    InsertAmplitude(Complex amp, size_t index);
    void    InsertAmplitude(std::vector<Complex> amp, size_t index);
    
    std::vector<Complex>
            GetAmplitude(const std::vector<int>& spins) const;
    std::vector<Complex>
            GetAmplitude(size_t index) const;
    
    double  SumSquare() const;
    
    void    Recreate(Amplitude_Tensor* newamps);
    void    CreateTrivial();
    
    size_t  Size() const;
    size_t  ColorSize() const;
    const CMatrix* GetColorMatrix() const;
    void    SetColorMatrix(const CMatrix* colormatrix);

    bool    Contains(const Particle* part) const;
    
    friend std::ostream& operator<<( std::ostream&, const Amplitude_Tensor &);
  };
}
#endif
