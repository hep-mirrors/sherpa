#ifndef PHASIC_Selectors_Selector_H
#define PHASIC_Selectors_Selector_H

#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Selector_List.H"

#include <limits>

namespace ATOOLS { class Data_Reader; }

namespace PHASIC {

  class Process_Base;

  class Selector_Log {
  private:

    std::string m_name;
    long int    m_rejected, m_passed;

  public :

    inline Selector_Log(const std::string &name): 
      m_name(name), m_rejected(0), m_passed(0) {}

    void Output();

    inline long int Rejections() const { return m_rejected; }
    inline long int Passed() const     { return m_passed;   }

    inline void ChangeName(std::string newname) { m_name=newname; }

    inline int  Hit(bool hit) 
    { 
      if (hit) { ++m_rejected; return 1; } 
      else { ++m_passed; return 0;}
    }

  };// end of class Selector_Log

  struct Selector_Key: public std::vector<std::vector<std::string> > {
    Process_Base        *p_proc;
    ATOOLS::Data_Reader *p_read;

    bool m_del;

    std::string m_key;

    inline Selector_Key(Process_Base *const proc,
			ATOOLS::Data_Reader *const read,const bool del=false):
      p_proc(proc), p_read(read), m_del(del) {}

    ~Selector_Key();

    void ReadData(const std::string& m_path, const std::string& m_file);
    void SetData(const std::string &tag,const std::vector<std::string> &args);

  };// end of struct Selector_Key

  std::ostream & operator<<(std::ostream & s,const PHASIC::Selector_Key & sk);

  class Selector_Base {
  protected :

    std::string m_name;
    bool        m_on,m_isnlo;

    Selector_Log *m_sel_log;
    Process_Base *p_proc;
    std::vector<Selector_Base*> m_sels;

    size_t m_nin, m_nout, m_n;
    int m_pass;

    ATOOLS::Flavour *p_fl;
    ATOOLS::NLO_subevt *p_sub;

    double m_smin,m_smax;

  public :

    Selector_Base(const std::string &name,Process_Base *const proc=NULL);

    virtual ~Selector_Base();

    // member functions
    static void ShowSyntax(const int mode);

    virtual bool Trigger(const ATOOLS::Vec4D_Vector &p,
                         const ATOOLS::Flavour *fl=NULL, size_t n=0);
    virtual bool RSTrigger(ATOOLS::NLO_subevtlist *const subs);
    virtual bool Trigger(ATOOLS::Selector_List &sl) = 0;

    virtual void BuildCuts(Cut_Data *) = 0;
    virtual void AddOnshellCondition(std::string,double);
    virtual void Output();

    void ReadInSubSelectors(const Selector_Key &key,size_t idx);

    // inline functions
    inline void SetProcess(Process_Base *const proc)
    { p_proc=proc; }
    inline const PHASIC::Process_Base *const Process(Process_Base *const proc) const
    { return p_proc; }

    inline std::string Name() const { return m_name; }

    inline int Pass() const { return m_pass; }

  };// end of class Selector_Base

  typedef ATOOLS::Getter_Function
  <Selector_Base,Selector_Key> Selector_Getter;

}

#endif
