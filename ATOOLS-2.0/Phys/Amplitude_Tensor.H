#ifndef Amplitude_Tensor_h
#define Amplitude_Tensor_h

#include "MyComplex.H"
#include "Particle.H"
#include "Blob.H"
#include "Matrix.H"
#include "Message.H"
#include <algorithm>

namespace ATOOLS
{
  bool SortByFirst(const std::pair<int,int> p1, const std::pair<int,int> p2);

  template<class Value>
  class Spin_Structure : public std::vector<Value> {
  protected:
    std::vector<int>    m_spins;
    
    size_t GetNumber(const std::vector<int>& spins) const
    {
      if(spins.size()!=m_spins.size()) {
	msg.Error()<<METHOD<<" Error: wrong size of spin std::vector."<<std::endl;
      }
      int mult(1);
      size_t num(0);
      for(size_t i=0; i<spins.size(); i++) {
	num += mult * spins[i];
	mult *= m_spins[i];
      }
      if(num>this->size()) {
	msg.Error()<<METHOD<<" Error: tried to access value out of bounce. "
	  <<"num="<<num<<" > "<<this->size()<<std::endl;
      }
      return num;
    };
    
    size_t GetNumber(std::vector<std::pair<int,int> >& spins) const
    {
      sort(spins.begin(),spins.end(),SortByFirst);
      
      if(spins.size()!=m_spins.size()) {
	msg.Error()<<METHOD<<" Error: wrong size of spin std::vector."<<std::endl;
      }
      int mult(1);
      size_t num(0);
      for(size_t i=0; i<spins.size(); i++) {
	num += mult * spins[i].second;
	mult *= m_spins[i];
      }
      if(num>this->size()) {
	msg.Error()<<METHOD<<" Error: tried to access value out of bounce. "
	  <<"num="<<num<<" > "<<this->size()<<std::endl;
      }
      return num;
    };
    
    std::vector<int> GetSpinCombination(size_t number) const
    {
      std::vector<int> spins(m_spins.size());
      for(size_t i=0;i<m_spins.size();i++) {
	spins[i] = number%m_spins[i];
	number = (number-spins[i])/m_spins[i];
      }
      return spins;
    };

  public:
    Spin_Structure(Flavour* flavs, size_t size)
    {
      m_spins = std::vector<int>(size);

      size_t n=1;
      for(size_t i=0;i<size;i++) {
        if(flavs[i].IsVector() && flavs[i].PSMass()==0) m_spins[i] = 2;
        else m_spins[i] = flavs[i].IntSpin()+1;
	n*=m_spins[i];
      }
      this->resize(n);
    };
    
    Spin_Structure(const Particle_Vector& particles)
    {
      m_spins = std::vector<int>(particles.size());

      size_t n=1;
      for(size_t i=0;i<particles.size();i++) {
        if(particles[i]->Flav().IsVector() && particles[i]->Flav().PSMass()==0) m_spins[i] = 2;
        else m_spins[i] = particles[i]->Flav().IntSpin()+1;
	n*=m_spins[i];
      }
      this->resize(n);
    };
    
    ~Spin_Structure() {};

    void    Insert(Value value, std::vector<std::pair<int,int> >& spins)
    {
      (*this)[GetNumber(spins)]=value;
    };
    
    void    Insert(Value value, size_t index)
    {
      (*this)[index]=value;
    };

    Value   Get(const std::vector<int>& spins) const
    {
      return (*this)[GetNumber(spins)];
    };
    
    Value   Get(size_t index) const
    {
      return (*this)[index];
    };

    void    CreateTrivial(Value value)
    {
      size_t n = this->size();
      this->clear();
      this->resize(n,value);
    };

    std::vector<int> Spins()
    {
      return m_spins;
    }
  };

//   typedef Spin_Structure<Complex> Spin_Amplitudes;
  class Spin_Amplitudes : public Spin_Structure<Complex> {
  public:
    Spin_Amplitudes(Flavour* flavs, size_t size);
    Spin_Amplitudes(const Particle_Vector& particles);
    double  SumSquare() const;
  };

  class Amplitude_Tensor : public Spin_Structure<std::vector<Complex> > {
    Particle_Vector                     m_particles;
    const CMatrix                     * p_colormatrix;
    Amplitude_Tensor                  * p_contracted;
    
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
    
    double  SumSquare() const;

    size_t  ColorSize() const;
    const CMatrix* GetColorMatrix() const;
    void    SetColorMatrix(const CMatrix* colormatrix);

    bool    Contains(const Particle* part) const;
    const Particle_Vector& Particles() const;

    friend std::ostream& operator<<( std::ostream&, const Amplitude_Tensor &);
  };

  /*! 
    \class Spin_Structure
    \brief Storing objects by spin combination.
    
    This class provides methods to access and manipulate arbitrary objects, by
    specifying the corresponding spin combination. It inherits from STL vector
    which contains one 'object' for each spin combination.
  */

  /*!
    \var std::vector<int> Spin_Structure::m_spins
    \brief vector containing the number of spin combinations for each particle "node".
  */

  /*!
    \var std::vector<Value>
    \brief storage of objects
    
    Objects are ordered as in this example (where the node part1 has
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
    \fn size_t Spin_Structure::GetNumber(const vector<int>& spins) const
    \brief Determine number of given combination in the STL vector.

    Here, the spins have to be specified in the same order as in m_spins (i.e.
    as the flavours or particles in the constructor where specified).
  */

  /*!
    \fn size_t Spin_Structure::GetNumber(vector<pair<int,int> >& spins) const
    \brief Determine number of given combination in the STL vector.

    Here, the spins (second int in pair) can be in arbitrary order, as long as
    the first int in the pair specifies the number in m_spins it belongs to.
  */

  /*!
    \fn std::vector<int> Spin_Structure::GetSpinCombination(size_t number) const
    \brief Determine spin combination from number of combination in the vector.

    The other way around from the GetNumber methods. Useful to
    traverse through all spin combinations.
  */

  /*!
    \fn Spin_Structure::Spin_Structure(Flavour* flavs, size_t size)
    \brief Constructor
  */

  /*!
    \fn Spin_Structure::Spin_Structure(Particle_Vector particles)
    \brief Constructor

    Particles' flavours determine the number of spin combinations of each node.
  */

  /*!
    \fn Spin_Structure::~Spin_Structure()
    \brief Destructor

    Doesn't do anything, because there are no pointer members.
  */

  /*!
    \fn void Spin_Structure::Insert(Value value, vector<pair<int,int> >& spins)
    \brief Inserts value at the right position determined by GetNumber.
  */

  /*!
    \fn void Spin_Structure::Insert(Complex amp, size_t index)
    \brief Inserts value at the given position.
  */

  /*!
    \fn Value Spin_Structure::Get(const vector<int>& spins) const
    \brief Retrieves the value from the right position determined by GetNumber.
  */

  /*!
    \fn vector<Complex> Spin_Structure::Get(size_t index) const;
    \brief Retrieves the value from the given position.
  */

  /*!
    \fn void Spin_Structure::CreateTrivial(Value value)
    \brief Inserts the given object for all spin combinations.
  */


  /*! 
    \class Amplitude_Tensor
    \brief Storing amplitudes for all spin and color combinations of one process.
    
    This class inherits from Spin_Structure to provide methods for accessing
    and manipulating amplitudes, by specifying the corresponding spin
    combination. For each spin combination it contains a vector(!) of complex
    numbers, each representing one color combination.
    It also keeps track of the particles involved with these amplitudes, so that
    you can easily contract two amplitude tensors over the spin index of one
    particle.
    To allow for sum-squaring over those amplitudes with the correct color
    factors, it also contains a color matrix, which should have the same rank
    as the number of amplitudes per spin combination.
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
    \fn double Amplitude_Tensor::SumSquare() const
    \brief Calculates the spin sum of all squared amplitudes (or of their summed
           color contraction, if necessary)
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
