#ifndef PHASIC_Selectors_Selector_H
#define PHASIC_Selectors_Selector_H

#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"

#include <limits>

namespace ATOOLS { class Data_Reader; }

namespace PHASIC {

  class Process_Integrator;

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

    inline int  Hit() { ++m_rejected; return 0; }
    inline int  Hit(bool hit) 
    { 
      if (hit) { ++m_rejected; return 1; } 
      else { ++m_passed; return 0;}
    }

  };// end of class Selector_Log

  struct Selector_Key: public std::vector<std::vector<std::string> > {
    Process_Integrator  *p_proc;
    ATOOLS::Data_Reader *p_read;

    bool m_del;

    std::string m_key;

    inline Selector_Key(Process_Integrator *const proc,
			ATOOLS::Data_Reader *const read,const bool del=false):
      p_proc(proc), p_read(read), m_del(del) {}

    ~Selector_Key();

    void ReadData(const std::string& m_path, const std::string& m_file);
    void SetData(const std::string &tag,const std::vector<std::string> &args);

  };// end of struct Selector_Key

  class Selector_Base {
  protected :

    std::string m_name;

    Selector_Log       *m_sel_log;
    Process_Integrator *p_proc;

    int m_n, m_nin, m_nout;

    ATOOLS::Flavour *m_fl;

    double m_smin,m_smax;

  public :

    inline Selector_Base(const std::string &name): 
      m_name(name), m_sel_log(NULL), p_proc(NULL),
      m_n(0), m_nin(0), m_nout(0), m_fl(NULL),
      m_smin(0.), m_smax(std::numeric_limits<double>::max()) {}

    virtual ~Selector_Base();

    // member functions
    static void ShowSyntax(const int mode);

    virtual bool Trigger(const ATOOLS::Vec4D_Vector &p) = 0;
    virtual bool JetTrigger(const ATOOLS::Vec4D_Vector &,
			    ATOOLS::NLO_subevtlist *const subs);
    virtual bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p);

    virtual void BuildCuts(Cut_Data *);
    virtual void AddOnshellCondition(std::string,double);
    virtual void Output();
    
    virtual int IsConditional();

    // inline functions
    inline void SetProcess(Process_Integrator *const proc) 
    { p_proc=proc; }

    inline std::string Name() const { return m_name; }

  };// end of class Selector_Base

  typedef ATOOLS::Getter_Function
  <Selector_Base,Selector_Key> Selector_Getter;

}

#endif
