//Message header
#ifndef Message_h
#define Message_h

#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include "STL_Tools.H"

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

    std::set<std::string,String_Sort> m_contextevents;
    std::set<std::string,String_Sort> m_contextinfo;
    std::set<std::string,String_Sort> m_contexttracking;
    std::set<std::string,String_Sort> m_contextdebugging;

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
    void Init(const std::string& level,const std::string &logfile="");
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

    bool LevelIsEvents(const std::string& context) const;
    bool LevelIsInfo(const std::string& context) const;
    bool LevelIsTracking(const std::string& context) const;
    bool LevelIsDebugging(const std::string& context) const;

    inline const std::set<std::string,String_Sort>& ContextEvents() {
      return m_contextevents;
    }
    inline const std::set<std::string,String_Sort>& ContextInfo() {
      return m_contextinfo;
    }
    inline const std::set<std::string,String_Sort>& ContextTracking() {
      return m_contexttracking;
    }
    inline const std::set<std::string,String_Sort>& ContextDebugging() {
      return m_contextdebugging;
    }
  };// end of class Message

  extern Message *msg;

  class Indentation {
  private:

    size_t m_col;

  public:

    // constructor
    inline Indentation(const size_t col=2) { msg->m_indent+=m_col=col; }

    // destructor
    inline ~Indentation() { msg->m_indent-=m_col; }

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

#define msg_LevelIsEvents() \
  (ATOOLS::msg->LevelIsEvents() || \
   (!ATOOLS::msg->ContextEvents().empty() && ATOOLS::msg->LevelIsEvents(__PRETTY_FUNCTION__)))
#define msg_LevelIsInfo() \
  (ATOOLS::msg->LevelIsInfo() || \
   (!ATOOLS::msg->ContextInfo().empty() && ATOOLS::msg->LevelIsInfo(__PRETTY_FUNCTION__)))
#define msg_LevelIsTracking() \
  (ATOOLS::msg->LevelIsTracking() || \
   (!ATOOLS::msg->ContextTracking().empty() && ATOOLS::msg->LevelIsTracking(__PRETTY_FUNCTION__)))
#define msg_LevelIsDebugging() \
  (ATOOLS::msg->LevelIsDebugging() || \
   (!ATOOLS::msg->ContextDebugging().empty() && ATOOLS::msg->LevelIsDebugging(__PRETTY_FUNCTION__)))

#define msg_Error() ATOOLS::msg->Error()
#define msg_Out() ATOOLS::msg->Out()
#define msg_LogFile() ATOOLS::msg->LogFile()
#define msg_Events() if (msg_LevelIsEvents()) ATOOLS::msg->Out()
#define msg_Info() if (msg_LevelIsInfo()) ATOOLS::msg->Out()
#define msg_Tracking() if (msg_LevelIsTracking()) ATOOLS::msg->Out()
#define msg_Debugging() if (msg_LevelIsDebugging()) ATOOLS::msg->Out()

#define msg_Indent() ATOOLS::Indentation indent
#define msg_Indentation(COL) ATOOLS::Indentation indent(COL)

#define bm_back ATOOLS::bm::back
#define bm_cr ATOOLS::bm::cr
#define bm_bell ATOOLS::bm::bell

#define mm_up(LINES) ATOOLS::mm(LINES,mm::up)
#define mm_down(LINES) ATOOLS::mm(LINES,mm::down)
#define mm_left(COLUMNS) ATOOLS::mm(COLUMNS,mm::left)
#define mm_right(COLUMNS) ATOOLS::mm(COLUMNS,mm::right)

#define METHOD \
  ATOOLS::msg->ExtractMethodName(__PRETTY_FUNCTION__)

#define PRINT_METHOD \
  std::cout<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset<<std::endl

#define PRINT_INFO(INFO) \
  std::cout<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset \
           <<":("<<ATOOLS::om::green<<"\""<<INFO<<"\"" \
           <<ATOOLS::om::reset<<")"<<std::endl

#define VAR(ARG) #ARG<<"="<<ARG

#define PRINT_VAR(VAR) \
  msg_Out()<<ATOOLS::om::blue<<METHOD<<":"<<__LINE__<<" "   \
  <<#VAR<<ATOOLS::om::reset<<"="<<ATOOLS::om::green \
  <<VAR<<ATOOLS::om::reset<<std::endl

#define PRINT_FUNC(a1) \
  msg_Out()<<ATOOLS::om::red<<METHOD<<"("<<ATOOLS::om::green \
  <<a1<<ATOOLS::om::red<<")"<<ATOOLS::om::reset<<std::endl

#define DEBUG_INFO(INFO) \
  msg_Debugging()<<ATOOLS::om::blue<<METHOD<<":"<<__LINE__      \
  <<"("<<ATOOLS::om::green<<"\""<<INFO<<"\""                    \
  <<ATOOLS::om::reset<<")"<<std::endl

#define DEBUG_VAR(VAR) \
  msg_Debugging()<<ATOOLS::om::blue<<METHOD<<":"<<__LINE__<<" " \
  <<#VAR<<ATOOLS::om::reset<<"="<<ATOOLS::om::green \
  <<VAR<<ATOOLS::om::reset<<std::endl

#define DEBUG_FUNC(a1) \
  msg_Debugging()<<ATOOLS::om::red<<METHOD<<"("<<ATOOLS::om::green \
  <<a1<<ATOOLS::om::red<<")"<<ATOOLS::om::reset<<std::endl; \
  if (msg_LevelIsDebugging()) msg_Indent()

#endif

