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

    void    Contract( Particle* part1, Amplitude_Tensor* const amps, Particle* part2 );
    double  SoftContract( Particle* part1, Amplitude_Tensor* const amps, Particle* part2 );

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

  /*! 
    \class Amplitude_Tensor
    \brief Storing amplitudes for all spin and color combinations of one process.
    
    This class provides methods to access and manipulate amplitudes, by
    specifying the corresponding spin combination. For each spin combination it
    contains a vector of complex numbers, each representing one color
    combination.
    It also keeps track of the particles involved with this amplitude, so that
    you can easily contract two amplitude tensors over the spin index of one
    particle.
    To allow for sum-squaring over those amplitudes with the correct color
    factors, it also contains a color matrix, which should have the same rank
    as the number of amplitudes per spin combination.
  */

  /*!
    \var std::vector<std::vector<Complex> >  Amplitude_Tensor::m_amplitudes
    \brief storage of complex amplitudes
    
    Amplitudes are considered to be ordered as in this example (where part1 has
    2 spin combinations, part2 has 1, part3 has 3, part4 has 2):
        |  part1 |    part2 |    part3 |    part4 |
      0 |      0 |        0 |        0 |        0 |
      1 |      1 |        0 |        0 |        0 |
      2 |      0 |        0 |        1 |        0 |
      3 |      1 |        0 |        1 |        0 |
      4 |      0 |        0 |        2 |        0 |
      5 |      1 |        0 |        2 |        0 |
      6 |      0 |        0 |        0 |        1 |
      7 |      1 |        0 |        0 |        1 |
      8 |      0 |        0 |        1 |        1 |
      9 |      1 |        0 |        1 |        1 |
     10 |      0 |        0 |        2 |        1 |
     11 |      1 |        0 |        2 |        1 |
  */

  /*!
    \var Particle_Vector Amplitude_Tensor::m_particles;
    \brief particles which are involved with the amplitude
  */

  /*!
    \var const CMatrix* Amplitude_Tensor::p_colormatrix;
    \brief pointer to the matrix which contains the color factors for
           sum-squaring the amplitude
  */

  /*!
    \fn size_t Amplitude_Tensor::GetAmplitudeNumber(const vector<int>& spins) const
    \brief Determine number of amplitude (in m_amplitudes) of given combination

    Here, the spins have to be specified in the same order as m_particles.
  */

  /*!
    \fn size_t Amplitude_Tensor::GetAmplitudeNumber(vector<pair<int,int> >& spins) const
    \brief Determine number of amplitude (in m_amplitudes) of given combination

    Here, the spins (second int in pair) can be in arbitrary order, as long as
    the first int in the pair specifies the number in m_particles it belongs to.
  */

  /*!
    \fn std::vector<int> Amplitude_Tensor::GetSpinCombination(size_t ampnumber) const
    \brief Determine spin combination from number of amplitude (in m_amplitudes)

    The other way around from the GetAmplitudeNumber methods. Useful to
    traverse through all spin combinations.
  */

  /*!
    \fn Amplitude_Tensor::Amplitude_Tensor(Particle_Vector particles)
    \brief Constructor

    You have to specify the particle pointers, because that way we can later
    find the amplitude of a certain spin combination, and furthermore even
    contract over a particle pointer, which is present in two Amplitude_Tensor's
  */

  /*!
    \fn Amplitude_Tensor::~Amplitude_Tensor()
    \brief Destructor

    Doesn't do anything, because the only pointer member is the colormatrix,
    which is set from outside, and not new'ed within Amplitude_Tensor.
  */

  /*!
    \fn void Amplitude_Tensor::InsertAmplitude(Complex amp, vector<pair<int,int> >& spins)
    \brief Inserts amplitude at the right position determined by
           GetAmplitudeNumber.
  */

  /*!
    \fn void Amplitude_Tensor::InsertAmplitude(Complex amp, vector<pair<int,int> >& spins)
    \brief Inserts an amplitude vector of different color combinations at the
           right position determined by GetAmplitudeNumber.
  */

  /*!
    \fn void Amplitude_Tensor::InsertAmplitude(Complex amp, size_t index)
    \brief Inserts a single amplitude at the given position.
  */

  /*!
    \fn void Amplitude_Tensor::InsertAmplitude(vector<Complex> amp, size_t index)
    \brief Inserts an amplitude vector of different color combinations at the
           given position.
  */

  /*!
    \fn vector<Complex> Amplitude_Tensor::GetAmplitude(const vector<int>& spins) const
    \brief Retrieves the amplitude vector of different color combinations from
           the right position determined by GetAmplitudeNumber.
  */

  /*!
    \fn vector<Complex> Amplitude_Tensor::GetAmplitude(size_t index) const;
    \brief Retrieves the amplitude vector of different color combinations from
           the given position.
  */

  /*!
    \fn double Amplitude_Tensor::SumSquare() const
    \brief Calculates the spin sum of all squared amplitudes (or of their summed
           color contraction, if necessary)
  */

  /*!
    \fn void Amplitude_Tensor::Recreate(Amplitude_Tensor* newamps)
    \brief Just a simple way to copy another Amplitude_Tensor into this one.
  */

  /*!
    \fn void Amplitude_Tensor::CreateTrivial()
    \brief Inserts a (1.0,0.0) for all amplitudes.
  */

  /*!
    \fn bool Amplitude_Tensor::Contains(const Particle* part) const
    \brief Method to check if a particle appears in this Amplitude_Tensor.
  */

  /*!
    \fn friend Amplitude_Tensor Contraction(
                                     Particle* part1, Particle* part2,
                                     Amplitude_Tensor* const amps1,
                                     Amplitude_Tensor* const amps2
                                   );
    \brief Contracts amps1 and amps2 over the spin indices of part1 and part2.

    part1 and part2 may be the same, but due to particles sometimes just passing
    through a blob (like Shower or Fragmentation of tau), they don't have to.
    Of course, in a physical meaning, they have to represent the same particle.
  */



}


#endif
