#ifndef Random_H
#define Random_H

#include <fstream>
#include "MathTools.H"

namespace AMATOOLS {
  class Random {
    long     m_id;
    int      m_inext,m_inextp;
    long     m_ma[56];
    
    // for read and write status
    long int m_written;
    char     m_outname[150];
    std::fstream * m_outstream;

    double   Ran1(long *idum);
    double   Ran2(long *idum);
    double   Ran3();
    void     InitRan3(long*);

  public:
    Random(long nid) {SetSeed(nid);}
    
    inline double Get() {return Ran3();}
    inline double GetNZ() {      // not zero random number
      double ran1;
      do {
	ran1 = Get();
      } 
      while (ran1==0.); 
      return ran1;
    }
    inline long GetSeed() {return m_id;}
    inline void SetSeed(long nid) {
      m_id = nid<0 ? nid : -nid;
      InitRan3(&m_id);
      // for status recovery
      m_written=0;    
      m_outstream=0;
      strcpy(m_outname,"");
    }
    
    inline double Theta() {
      return acos(2.*Get()-1.);
    }  

    int WriteOutStatus(const char*);

    void ReadInStatus(const char* filename, long int index=0 );

    ~Random() { 
      if (m_outstream!=0) {
	m_outstream->close();
	delete m_outstream;
      }
    } 
  };

  extern Random ran;


  // --------------------------------------------------
  //         Doxygen part
  // --------------------------------------------------

  /*!
    \file
    \brief contains the class AMATOOLS::Random
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

  

}


#endif
