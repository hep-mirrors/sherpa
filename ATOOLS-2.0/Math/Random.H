#ifndef Random_H
#define Random_H

#include <fstream>
#include "MathTools.H"
#include "Exception.H"

namespace ATOOLS {

  class Random: public Terminator_Object {
  private:

    long int m_id, m_ma[56], m_written;
    int      m_inext,m_inextp;
    char     m_outname[150];

    std::fstream *p_outstream;

    double Ran1(long *idum);
    double Ran2(long *idum);
    double Ran3();

    void InitRan3(long*);

    void PrepareTerminate();

  public:

    // constructor
    Random(long nid): p_outstream(NULL) { SetSeed(nid); }

    // destructor
    ~Random();
    
    // member functions
    void SetSeed(long nid);

    inline long int GetSeed() { return m_id; }

    int  WriteOutStatus(const char *outfile);
    void ReadInStatus(const char *infile,long int index=0);

    inline double Get()   { return Ran2(&m_id);       }
    inline double Theta() { return acos(2.*Get()-1.); }  
    double        GetNZ();

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
    \fn void   Random::InitRan3(long*)
    \brief initializes the random number generator
  */

  /*!
    \fn Random::Random(long nid) 
    \brief Constructor initialises the random number generator with a given seed
  */

  /*!
    \fn Random::~Random() 
    \brief Destructor
  */

  /*!
    \fn inline double Random::Get() 
    \brief is the main routine, returns a single random number using Ran3()
  */

  /*!
    \fn inline double Random::GetNZ()
    \brief retrun a not zero random number
  */

  /*!
    \fn inline long Random::GetSeed()
    \brief returns a the seed
  */

  /*!
    \fn inline void Random::SetSeed(long nid)
    \brief sets a new seed and (re)initializes the random number generator
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
