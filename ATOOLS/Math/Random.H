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

  class Random: public Terminator_Object {
  private:

    int      activeGenerator;

    long int m_id, m_ma[56];
    long int m_sid, m_sma[56];
    int      m_inext,m_inextp;
    int      m_sinext,m_sinextp;

    External_RNG *p_external;

    // Variables for Ran4()
    struct Ran4Status{
      char idTag[16];
      int i97,j97;
      double u[97],c,cd,cm; } status;
    Ran4Status backupStat;   // Backup data for Save-/RestoreStatus()

    double Ran1(long *idum);
    double Ran2(long *idum);
    double Ran3();
    double Ran4();              

    void InitRan3(long*);
 
    bool ReadInStatus(const std::string &path);
    void PrepareTerminate();

    // temporary methods for Ran4()
    int  WriteOutStatus4(const char *outfile); 
    void ReadInStatus4(const char * filename);
    void SaveStatus4();
    void RestoreStatus4();

  public:

    // constructors
    Random(long nid);  // initialization for Ran2()
    Random(int, int);  // initialization for Ran4()

    // destructor
    ~Random();
    
    // member functions
    bool InitExternal(const std::string &path,const std::string &file);
    void SetSeed(long nid);  // seed for Rnd2()
    void SetSeed(int, int);  // seed for Rnd4()
    inline long int GetSeed() { return m_id; }

    int  WriteOutStatus(const char *outfile);
    int  WriteOutStatus(std::ostream &os,const size_t &idx);
    void ReadInStatus(const char *infile);
    size_t ReadInStatus(std::istream &is,const size_t &idx);
    void SaveStatus();
    void RestoreStatus();

    // return uniformly distributed random number in [0,1] using active Generator
    double Get();   
    ptrdiff_t operator() (ptrdiff_t max);
    // produce two Gaussian distributed random number using active Generator
    void   Gaussian(double &,double &);

    inline double Theta() { return acos(2.*Get()-1.); }  
    double        GetNZ();
    
    inline External_RNG* GetExternalRng() { return p_external; }
    

  };// end of class Random

  extern Random ran;

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
    \fn double Random::Ran1(long *idum)
    \brief is a very fast but simple random number routine
  */

  /*!
    \fn double Random::Ran2(long *idum)
    \brief is a very stable and powerful random number routine
  */

  /*!
    \fn double Random::Ran3()
    \brief is a good and fast random number routine
  */

  /*!
    \fn double Random::Ran4()
    \brief a new random generator that still needs to be tested
  */

  /*!
    \fn void   Random::InitRan3(long*)
    \brief initializes the random number generator
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
