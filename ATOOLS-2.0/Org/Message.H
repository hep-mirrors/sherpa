//Message header
#ifndef Message_h
#define Message_h

#include <iostream>
#include <fstream>

namespace ATOOLS {

  struct om {

    enum code {
      none     =      0,
      reset    =      1,
      bold     =      2,
      blink    =      4,
      underln  =      8,
      blackbg  =     16,
      red      =     32,
      green    =     64,
      blue     =    128,
      brown    =    256,
      violet   =    512,
      lblue    =   1024,
      grey     =   2048,  
      redbg    =   4096,
      greenbg  =   8192,
      bluebg   =  16384,
      brownbg  =  32768,
      violetbg =  65536,
      lbluebg  = 131072,
      greybg   = 262144  
    };

  };// end of struct om

  std::ostream &operator<<(std::ostream &str,const om::code modifier);

  class Message {
  private:

    std::ofstream *m_no; 
    std::ostream  *m_output, *m_error;

    int m_file, m_level;

    bool m_modifiable;

  public:

    // constructor
    Message();

    // member functions
    void Init(const int level);
    void SetStandard();
    void SetFile(const char *file);
    void SetPrecision(const int precision);

    // inline functions
    inline void SetLevel(const int level)     { m_level = level;    }
    inline void SetModifiable(const bool mod) { m_modifiable = mod; }

    inline int  Level()      { return m_level;      }
    inline bool Modifiable() { return m_modifiable; }

    inline void SetNoStream(std::ostream *no)      { m_no = (std::ofstream*)no; }
    inline void SetOutStream(std::ostream *output) { m_output = output;         }
    inline void SetErrStream(std::ostream *error)  { m_error = error;           }

    inline std::ostream &Out()       { return *m_output; }

    inline std::ostream &Error()     { if (m_level >= 0) return *m_output; return *m_no; }
    inline std::ostream &Events()    { if (m_level & 1)  return *m_output; return *m_no; }
    inline std::ostream &Info()      { if (m_level & 2)  return *m_output; return *m_no; }
    inline std::ostream &Tracking()  { if (m_level & 4)  return *m_output; return *m_no; }
    inline std::ostream &Debugging() { if (m_level & 8)  return *m_output; return *m_no; }

  };// end of class Message

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
     - Info()  less important run information 
          (e.g. calculated partial cross sections) 
     - Tracking() detailed information of the integration process
     - Events()   output of each single event (event generation mode only)
     - Error()    critical error messages (they usually lead to a stop of 
           the program)
    .
  */
}

#endif

