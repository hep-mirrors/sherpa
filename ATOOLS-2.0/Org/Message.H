//Message header
#ifndef Message_h
#define Message_h

#include <iostream>
#include <fstream>

namespace ATOOLS {

  namespace om {

    extern const char *reset;
    extern const char *bold;
    extern const char *blink;
    extern const char *backgd;
    extern const char *red;
    extern const char *green;
    extern const char *blue;
    extern const char *brown;
    extern const char *violet;
    extern const char *lblue;
    extern const char *grey;
    
  }// end of namespace om

  class Message {
    std::ofstream * no; 
    std::ostream  * output;
    std::ostream  * error;
    int file;
    int level;
  public:
    //Standard-output
    Message();
    void Init(int _level);

    int  Level() { return level; }

    void SetLevel(int _level) { level = _level; }
    void SetFile(char *);
    void SetStandard();
    void SetNoStream(std::ostream *);
    void SetOutStream(std::ostream *);
    void SetErrStream(std::ostream *);
    void SetPrecision(int);

    std::ostream & Out() { return * output; }
    std::ostream & Error() { 
      if (level >= 0) return * output;
      return * no;
    }
    std::ostream & Events() { 
      if (level >= 1)  return * output;
      return * no;
    }
    std::ostream & Tracking() { 
      if (level >= 2)  return * output;
      return * no;
    }
    std::ostream & Debugging()  { 
      if (level >= 3)  return * output;
      return * no;
    }
  };

  extern Message msg;

  /*!
    \file
    \brief contains the class Message
  */

  /*!
    \class Message
    \brief This is the main output class.

    It is defined to easy manipulate the amount of produced output.
    There are four pipes that can be used.
     - Out()      Important run information 
          (e.g. start/end of an integration) and warnings 
     - Results()  less important run information 
          (e.g. calculated partial cross sections) 
     - Events() detailed information of the integration process
              and an output of each single event (event generation mode only)
     - Error()   critical error messages (they usually lead to a stop of 
           the program)
    .
  */
}

#endif

