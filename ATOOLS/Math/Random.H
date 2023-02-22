#ifndef ATOOLS_Math_Random_H
#define ATOOLS_Math_Random_H

#include <limits>
#include <fstream>
#include <stddef.h>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  struct RNG_Key {
  };

  class External_RNG {
  public:

    virtual ~External_RNG();

    virtual double Get() = 0;

    // Derived classes can opt in to support status saving/restoring.
    // In that case, override all three of the below functions and
    // in particular return true in CanRestoreStatus().
    virtual bool CanRestoreStatus() const { return false; }
    virtual void SaveStatus() {}
    virtual void RestoreStatus() {}

  };

  typedef Getter_Function<External_RNG,RNG_Key> RNG_Getter;

  class Marsaglia;
  class Sobol;

  class Random: public Terminator_Object {
  private:

    int      activeGenerator, m_npe;

    long int m_id, m_sid;

    std::stringstream m_lastincrementedseed;
    size_t m_nsinceinit, m_increment;

    External_RNG *p_external;

    Marsaglia *p_ran4[2];

    std::vector<double> m_p;

    double Ran2(long *idum);
    double Ran4();

    bool ReadInStatus(const std::string &path);
    void PrepareTerminate();

    // temporary methods for Ran4()
    int  WriteOutStatus4(const char *outfile);
    int  WriteOutStatus4(std::ostream &os,const size_t &idx);
    int  WriteOutSavedStatus4(const char *outfile);
    void ReadInStatus4(const char * filename);
    size_t ReadInStatus4(std::istream &is,const size_t &idx);
    void SaveStatus4();
    void RestoreStatus4();

  public:

    // constructors
    Random(long nid);  // initialization for Ran2()
    Random(unsigned int i1,unsigned int i2,unsigned int i3,
	   unsigned int i4,unsigned int i5,unsigned int i6);

    // destructor
    ~Random();

    // member functions
    bool InitExternal();
    void SetSeed(long nid);  // seed for Rnd2()
    void SetSeed(unsigned int i1,unsigned int i2,
		 unsigned int i3,unsigned int i4);
    long int GetSeed() { return m_id; }

    int  WriteOutStatus(const char *outfile);
    int  WriteOutStatus(std::ostream &os,const size_t &idx);
    int  WriteOutSavedStatus(const char *outfile);
    void ReadInStatus(const char *infile);
    size_t ReadInStatus(std::istream &is,const size_t &idx);
    bool CanRestoreStatus() const;
    void SaveStatus();
    void RestoreStatus();
    void FastForward(const size_t &n);
    void ResetToLastIncrementedSeed();
    void EraseLastIncrementedSeed()
    { m_lastincrementedseed.str(std::string()); }
    void SetSeedStorageIncrement(size_t inc) { m_increment=inc; }

    // return uniformly distributed random number in [0,1] using active Generator
    double Get();
    // produce Gaussian distributed random number using active Generator
    // according to the Marsaglia method
    double GetGaussian() {
      static auto hasSpare = false;
      static double spare;
      if (hasSpare) {
        hasSpare = false;
        return spare;
      }
      double ran1, ran2, R;
      do {
        ran1 = 2.*Get()-1.;
        ran2 = 2.*Get()-1.;
        R    = ran1*ran1+ran2*ran2;
      } while (R>1. || R==0.);
      R  = sqrt(-2.*log(R)/R);
      hasSpare = true;
      spare = ran2 * R;
      return ran1 * R;
    }
    // produce Poissonian distributed random number using active Generator
    double Poissonian(const double & lambda) {
      if(lambda>500.) {
	double u = Get();
	double v = Get();
	return int(sqrt(lambda)*sqrt(-2.*log(u))*cos(2.*M_PI*v)+lambda);
      }
      double disc(std::exp(-lambda)),p(1.);
      int N(0);
      while ((p*=Get())>disc) N++;
      return N;
    }
    double Theta() { return acos(2.*Get()-1.); }
    double GetNZ();

    External_RNG* GetExternalRng() { return p_external; }

    inline void StartCount() { m_npe=0; }
    inline void StopCount() { m_npe=-m_npe; }
    inline int  Count() { return m_npe; }

    inline void Clear() { m_p.clear(); }
    // adopt requirements for UniformRandomBitGenerator
    // (e.g. for use in std::shuffle)
    typedef size_t result_type;
    result_type operator()();
    constexpr static result_type min() noexcept {
      return std::numeric_limits<result_type>::min();
    }
    constexpr static result_type max() noexcept {
      return std::numeric_limits<result_type>::max();
    }

  };// end of class Random

  extern Random *ran;

  // --------------------------------------------------
  //         Doxygen part
  // --------------------------------------------------

  /*!
    \file
    \brief contains the class Random
  */

  /*!
    \class Random
    \brief supplies uniformly distributed random numbers
  */

  /*!
    \fn double Random::Ran2(long *idum)
    \brief is a very stable and powerful random number routine
  */

  /*!
    \fn double Random::Ran4()
    \brief a new random generator that still needs to be tested
  */

  /*!
    \fn Random::Random(long nid)
    \brief Constructor initialises the random number generator with a given seed
  */

  /*!
    \fn Random::Random(int, int)
    \brief A constructor that initializes the Rnd4() routine.

    Even though there are two different constructors for Rnd2() and Rnd4(),
    it is possible to switch between the two routines by calling the
    corresponding SetSeed() method.
  */

  /*!
    \fn Random::~Random()
    \brief Destructor
  */

  /*!
    \fn double Random::Get()
    \brief is the main routine, returns a single random number in [0,1]

    The number is determined either by using Ran2() or Ran4(), depending on
    which of the two generators is set in the activeGenerator variable.
  */

  /*!
    \fn double Random::GetNZ()
    \brief retrun a not zero random number
  */

  /*!
    \fn long Random::GetSeed()
    \brief returns a the seed

    No corresponding method for Ran4() exists so far.
  */

  /*!
    \fn void Random::SetSeed(long nid)
    \brief sets a new seed and (re)initializes the random number generator Rnd2()
  */

  /*!
    \fn void Random::SetSeed(int, int)
    \brief sets a new seed and (re)initializes the random number generator Rnd4()
  */

  /*!
    \fn double Random::Theta()
    \brief returns an angle \f$\phi\f$ for a uniform \f$cos(\phi)\f$ distribution
  */

  /*!
    \fn int Random::WriteOutStatus(char* filename)
    \brief writes the complete status the random generator in a file

    This method can be used to save the status of a random generator in a file
    the number of its entry in this file is return and can be used to read in
    the status via ReadInStatus().
  */

  /*!
    \fn void Random::ReadInStatus(char* filename, long int index=0 )
    \brief reads in a status from a file
  */

}// end of namespace ATOOLS

#endif
