//Message header
#ifndef Message_h
#define Message_h

#include <iostream>
#include <fstream>

namespace AORGTOOLS {
  class Message {
    std::ofstream *no; 
    std::ostream * output;
    std::ostream * error;
    int file;
    int level;
  public:
    //Standard-output
    Message();
    void Init(int _level);

    void SetLevel(int _level) { level = _level; }
    void SetFile(char *);
    void SetStandard();
    void SetNoStream(std::ostream *);
    void SetOutStream(std::ostream *);
    void SetErrStream(std::ostream *);
    void SetPrecision(int);

    std::ostream & Out();
    std::ostream & Error();
    std::ostream & Events();
    std::ostream & Tracking();
    std::ostream & Debugging();
  };

  /*
  inline std::ostream& Message::Out() { std::cout<<"o."<<endl;return std::cout; }
  
  inline std::ostream& Message::Error()  { std::cout<<"err."<<endl;return std::cout; }
  
  inline std::ostream& Message::Events() { std::cout<<"ev."<<endl;return std::cout; }

  inline std::ostream& Message::Tracking()  { std::cout<<"tr."<<endl;return std::cout; }
  
  inline std::ostream& Message::Debugging() { std::cout<<"dg."<<endl;return std::cout; }
  */


  inline std::ostream & Message::Out() { 
    return * output;
  }
  
  inline std::ostream & Message::Error() { 
    if (level >= 0) return * output;
    return * no;
  }
  
  inline std::ostream & Message::Events() { 
    if (level >= 1)  return * output;
    return * no;
  }
  
  inline std::ostream & Message::Tracking() { 
    if (level >= 2)  return * output;
    return * no;
  }
  
  inline std::ostream & Message::Debugging() { 
    if (level >= 3)  return * output;
    return * no;
  }

  extern Message msg;

  /*!
    \file
    \brief contains the class AORGTOOLS::Message
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

