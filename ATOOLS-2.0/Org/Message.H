//Message header
#ifndef Message_h
#define Message_h

#include <string>
#include <iostream>
#include <fstream>

namespace ATOOLS {

  struct bm {

    enum code {
      none = 0,
      back = 1,
      cr   = 2,
      bell = 4
    };

  };// end of struct bm

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

  struct mm {

    enum code {
      none  = 0,
      up    = 1,
      down  = 2,
      left  = 4,
      right = 8
    };

    int m_num, m_code;

    // constructor
    inline mm(const int num,const code type): 
      m_num(num), m_code(type) {}

  };// end of struct bm

  struct tm {

    enum code {
      none   = 0,
      curon  = 1,
      curoff = 2
    };

  };// end of struct tm

  std::ostream &operator<<(std::ostream &str,const bm::code modifier); 
  std::ostream &operator<<(std::ostream &str,const om::code modifier);

  std::ostream &operator<<(std::ostream &str,const mm modifier); 

  std::ostream &operator<<(std::ostream &str,const tm::code modifier); 

  class Indentation;

  class Message {
  private:

    std::ofstream *p_no, *p_logfile; 
    std::ostream  *p_output, *p_error;

    int m_file, m_level;

    int m_modifiable;

    size_t m_indent;

    void InitLogFile(const std::string &logfile); 

    inline const std::string Indent() const 
    { return std::string(m_indent,' '); }

    friend class Indentation;

  public:

    // constructor
    Message();

    // destructor
    ~Message(); 

    // member functions
    void Init(const int level,const std::string &logfile="");
    void SetPrecision(const int precision);
    void SetStandard();

    std::ostream &Out() const;

    std::ostream &Error() const;     
    std::ostream &Events() const;    
    std::ostream &Info() const;      
    std::ostream &Tracking() const;  
    std::ostream &Debugging() const; 

    std::ostream &LogFile() const;   

    std::string ExtractMethodName(std::string cmethod) const;

    // inline functions
    inline void SetLevel(const int level)     { m_level = level;    }
    inline void SetModifiable(const bool mod) { m_modifiable = mod; }

    inline int  Level() const      { return m_level;      }
    inline bool Modifiable() const { return m_modifiable; }

    inline void SetNoStream(std::ostream *const no)      { p_no=(std::ofstream*)no; }
    inline void SetOutStream(std::ostream *const output) { p_output=output;         }
    inline void SetErrStream(std::ostream *const error)  { p_error=error;           }

    inline bool LevelIsError() const     { return (m_level >= 0); }
    inline bool LevelIsEvents() const    { return (m_level & 1);  }
    inline bool LevelIsInfo() const      { return (m_level & 2);  }
    inline bool LevelIsTracking() const  { return (m_level & 4);  }
    inline bool LevelIsDebugging() const { return (m_level & 8);  }

  };// end of class Message

  extern Message msg;

  class Indentation {
  private:

    size_t m_col;

  public:

    // constructor
    inline Indentation(const size_t col=2) { msg.m_indent+=m_col=col; }

    // destructor
    inline ~Indentation() { msg.m_indent-=m_col; }

  };//end of class Indentation

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

#define msg_Events() if (ATOOLS::msg.LevelIsEvents()) ATOOLS::msg.Out()
#define msg_Info() if (ATOOLS::msg.LevelIsInfo()) ATOOLS::msg.Out()
#define msg_Tracking() if (ATOOLS::msg.LevelIsTracking()) ATOOLS::msg.Out()
#define msg_Debugging() if (ATOOLS::msg.LevelIsDebugging()) ATOOLS::msg.Out()

#define msg_Indent() ATOOLS::Indentation indent
#define msg_Indentation(COL) ATOOLS::Indentation indent(COL)

#define METHOD \
  ATOOLS::msg.ExtractMethodName(__PRETTY_FUNCTION__)

#define PRINT_METHOD \
  std::cout<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset<<std::endl

#define PRINT_INFO(INFO) \
  std::cout<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset \
           <<":("<<ATOOLS::om::green<<"\""<<INFO<<"\"" \
           <<ATOOLS::om::reset<<")"<<std::endl

#endif

