#ifndef ATOOLS_Org_Message_H
#define ATOOLS_Org_Message_H
//Message header

#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/CXXFLAGS.H"

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

  struct cm {

    enum code {
      none   = 0,
      curon  = 1,
      curoff = 2
    };

  };// end of struct cm

  std::ostream &operator<<(std::ostream &str,const bm::code modifier); 
  std::ostream &operator<<(std::ostream &str,const om::code modifier);

  std::ostream &operator<<(std::ostream &str,const mm modifier); 

  std::ostream &operator<<(std::ostream &str,const cm::code modifier); 

  class Indentation;

  class indentbuf : public std::streambuf {
  public:
    indentbuf(std::streambuf* basebuf);
    ~indentbuf();
    void Indent(size_t i=2);
    void DeIndent(size_t i=2);
    inline std::streambuf* BaseBuf() { return m_basebuf; }
    inline void SetBaseBuf(std::streambuf *buf) { m_basebuf=buf; }

  protected:
    std::streambuf* m_basebuf;
    size_t m_indent;
    bool at_start;

    int_type overflow(int_type ch);
    inline int sync() { return m_basebuf->pubsync(); }
  };

  class Message {
  private:

    std::ofstream m_devnull;
    std::ofstream *p_log;
    std::ostream m_output, m_error;
    indentbuf m_buf;

    int m_level;
    std::string m_logfile;

    std::set<std::string,String_Sort> m_contextevents;
    std::set<std::string,String_Sort> m_contextinfo;
    std::set<std::string,String_Sort> m_contexttracking;
    std::set<std::string,String_Sort> m_contextdebugging;
    std::set<std::string,String_Sort> m_contextiodebugging;

    int m_modifiable, m_mpimode;

  public:

    // constructor
    Message();

    // destructor
    ~Message(); 

    // member functions
    void Init();
    void SetStandard();

    std::ostream &Out();
    std::ostream &Error();
    std::ostream &Events();
    std::ostream &Info();
    std::ostream &Tracking();
    std::ostream &Debugging();
    std::ostream &IODebugging();

    std::string ExtractMethodName(std::string cmethod) const;

    // inline functions
    inline void SetLevel(const int level)     { m_level = level;    }
    inline void SetModifiable(const bool mod) { m_modifiable = mod; }
    inline void SetPrecision(const int precision) {
      m_output.precision(precision);
    }

    inline int  Level() const      { return m_level;      }
    inline bool Modifiable() const { return m_modifiable; }
    inline size_t Precision() const {
      return m_output.precision();
    }


    inline void Indent(size_t col) { m_buf.Indent(col); }
    inline void DeIndent(size_t col) { m_buf.DeIndent(col); }

    inline void SetOutStream(std::ostream& output) {
      m_output.rdbuf(output.rdbuf());
    }
    inline void SetOutStream(ATOOLS::indentbuf& buf) {
      m_output.rdbuf(&buf);
    }
    inline void SetErrStream(std::ostream& error) {
      m_error.rdbuf(error.rdbuf());
    }
    inline void SetErrStream(ATOOLS::indentbuf& buf) {
      m_error.rdbuf(&buf);
    }

    inline const std::string& LogFile() const { return m_logfile; }

    inline bool LevelIsError() const     { return (m_level >= 0); }
    inline bool LevelIsEvents() const    { return (m_level & 1);  }
    inline bool LevelIsInfo() const      { return (m_level & 2);  }
    inline bool LevelIsTracking() const  { return (m_level & 4);  }
    inline bool LevelIsDebugging() const { return (m_level & 8);  }
    inline bool LevelIsIODebugging() const { return (m_level & 32);  }

    bool LevelIsEvents(const std::string& context) const;
    bool LevelIsInfo(const std::string& context) const;
    bool LevelIsTracking(const std::string& context) const;
    bool LevelIsDebugging(const std::string& context) const;
    bool LevelIsIODebugging(const std::string& context) const;

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
    inline const std::set<std::string,String_Sort>& ContextIODebugging() {
      return m_contextiodebugging;
    }
  };// end of class Message

  extern Message *msg;

  class Indentation {
  private:

    size_t m_col;
    int m_mode;

  public:

    // constructor
    inline Indentation(const size_t col=2) : m_col(col), m_mode(0) {
    }

    inline void Activate(int mode) {
      m_mode=mode;
      if (m_mode & 2) msg->Out()<<om::red<<"{"<<om::reset<<std::endl;
      if (m_mode & 1) msg->Indent(m_col);
    }


    // destructor
    inline ~Indentation() {
      if (m_mode & 1) msg->DeIndent(m_col);
      if (m_mode & 2) msg->Out()<<om::red<<"}"<<om::reset<<std::endl;
    }

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
#define msg_LevelIsIODebugging() \
  (ATOOLS::msg->LevelIsIODebugging() || \
   (!ATOOLS::msg->ContextIODebugging().empty() && ATOOLS::msg->LevelIsIODebugging(__PRETTY_FUNCTION__)))

#define msg_Error() ATOOLS::msg->Error()
#define msg_Out() ATOOLS::msg->Out()
#define msg_Events() if (msg_LevelIsEvents()) ATOOLS::msg->Out()
#define msg_Info() if (msg_LevelIsInfo()) ATOOLS::msg->Out()
#define msg_Tracking() if (msg_LevelIsTracking()) ATOOLS::msg->Out()
#define msg_Debugging() if (msg_LevelIsDebugging()) ATOOLS::msg->Out()
#define msg_IODebugging() if (msg_LevelIsIODebugging()) ATOOLS::msg->Out()

#define msg_Indent() ATOOLS::Indentation indent; indent.Activate(1)
#define msg_Indentation(COL) ATOOLS::Indentation indent(COL);indent.Activate(1);

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
  msg_Out()<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset<<std::endl

#define PRINT_INFO(INFO) \
  msg_Out()<<ATOOLS::om::blue<<METHOD<<ATOOLS::om::reset        \
           <<":("<<ATOOLS::om::green<<"\""<<INFO<<"\"" \
           <<ATOOLS::om::reset<<")"<<std::endl

#define VAR(ARG) #ARG<<"="<<ARG

#define PRINT_VAR(VAR) \
  msg_Out()<<ATOOLS::om::blue<<__LINE__<<":"   \
  <<#VAR<<ATOOLS::om::reset<<"="<<ATOOLS::om::green \
  <<VAR<<ATOOLS::om::reset<<std::endl

#define PRINT_FUNC(a1) \
  msg_Out()<<ATOOLS::om::red<<METHOD<<"("<<ATOOLS::om::green \
  <<a1<<ATOOLS::om::red<<")"<<ATOOLS::om::reset; \
  ATOOLS::Indentation indent; \
  indent.Activate(3)

#define DEBUG_INFO(INFO) \
  msg_Debugging()<<ATOOLS::om::blue<<METHOD<<":"<<__LINE__      \
  <<"("<<ATOOLS::om::green<<"\""<<INFO<<"\""                    \
  <<ATOOLS::om::reset<<")"<<std::endl

#define DEBUG_VAR(VAR) \
  msg_Debugging()<<ATOOLS::om::blue<<__LINE__<<":" \
  <<#VAR<<ATOOLS::om::reset<<"="<<ATOOLS::om::green \
  <<VAR<<ATOOLS::om::reset<<std::endl

#define DEBUG_FUNC(a1) \
  msg_Debugging()<<ATOOLS::om::red<<METHOD<<"("<<ATOOLS::om::green \
  <<a1<<ATOOLS::om::red<<") "<<ATOOLS::om::reset; \
  ATOOLS::Indentation indent; \
  if (msg_LevelIsDebugging()) indent.Activate(3)

#endif

