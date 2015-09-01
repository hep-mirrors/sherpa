#ifndef ATOOLS_Math_Random_H
#define ATOOLS_Math_Random_H

#include <fstream>
#include <stddef.h>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  struct RNG_Key {
  };

  class External_RNG {
  public:

    virtual ~External_RNG();

    virtual double Get() = 0;

  };

  typedef Getter_Function<External_RNG,RNG_Key> RNG_Getter;

  class Marsaglia;

  class Random: public Terminator_Object {
  private:

    int      activeGenerator;

    long int m_id, m_sid;

    std::stringstream m_lastincrementedseed;
    size_t m_nsinceinit, m_increment;

    External_RNG *p_external;

    Marsaglia *p_ran4[2];

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
    bool InitExternal(const std::string &path,const std::string &file);
    void SetSeed(long nid);  // seed for Rnd2()
    void SetSeed(unsigned int i1,unsigned int i2,
		 unsigned int i3,unsigned int i4);
    inline long int GetSeed() { return m_id; }

    int  WriteOutStatus(const char *outfile);
    int  WriteOutStatus(std::ostream &os,const size_t &idx);
    int  WriteOutSavedStatus(const char *outfile);
    void ReadInStatus(const char *infile);
    size_t ReadInStatus(std::istream &is,const size_t &idx);
    void SaveStatus();
    void RestoreStatus();
    void FastForward(const size_t &n);
    void ResetToLastIncrementedSeed();
    inline void EraseLastIncrementedSeed() 
    { m_lastincrementedseed.str(std::string()); }
    inline void SetSeedStorageIncrement(size_t inc) { m_increment=inc; }

    // return uniformly distributed random number in [0,1] using active Generator
    double Get();   
    ptrdiff_t operator() (ptrdiff_t max);
    // produce two Gaussian distributed random number using active Generator
    void   Gaussian(double &,double &);
    // produce Poissonian distributed random number using active Generator
    inline double Poissonian(const double & lambda) {
      if(lambda>500.) {
	double u = Get();
	double v = Get();
	return int(lambda*sqrt(-2.*log(u))*cos(2.*M_PI*v)+lambda);
      }
      double disc(std::exp(-lambda)),p(1.);
      int N(0);
      while ((p*=Get())>disc) N++; 
      return N;
    }

    inline double Theta() { return acos(2.*Get()-1.); }  
    double        GetNZ();
    
    inline External_RNG* GetExternalRng() { return p_external; }
    

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
    \fn inline double Random::Get() 
    \brief is the main routine, returns a single random number in [0,1]

    The number is determined either by using Ran2() or Ran4(), depending on
    which of the two generators is set in the activeGenerator variable.
  */

  /*!
    \fn inline double Random::GetNZ()
    \brief retrun a not zero random number
  */

  /*!
    \fn inline long Random::GetSeed()
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
    \fn inline double Random::Theta()
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
