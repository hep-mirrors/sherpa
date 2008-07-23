#ifndef Selector_Base_h
#define Selector_Base_h

#include "Vector.H"
#include "Flavour.H"
#include "Message.H"
#include "Cut_Data.H"

#include <vector>

namespace ATOOLS {
  struct Mom_Data {
    std::vector<Flavour> flavs;
    std::vector<std::pair<double, double> >  bounds;
    std::string helps;
    int    help; 
    int    type;
  };

  class Selector_Data {
  protected :
    std::vector<Mom_Data> data;
  private :
    bool ReadInData(std::string path, std::string filename);
  public :
    Selector_Data(std::string path, std::string filename);
    inline Selector_Data(Selector_Data & sel) : data(sel.data) {}
    Selector_Data();
    ~Selector_Data()     {};
    int  Size()          { return data.size(); }
    Mom_Data RemoveData(int &);
    inline Mom_Data GetData(int i) { return data[i]; }
    void SetData(int ,const std::vector<Flavour> &,
		 std::vector<std::pair<double,double> > &,int);
    void AddData(int ,const std::vector<Flavour> &,
		 std::vector<std::pair<double,double> > &,int);
    void FillData(int,int &,std::vector<Flavour> &,
		  std::vector<std::pair<double,double> > &,int &);
    void ControlOutput();
  };


  class Selector_Log;

  class Selector_Base {
  protected :
    int            m_n, m_nin, m_nout;
    Flavour      * m_fl;
    double         m_smin,m_smax;
    std::string    m_name, m_procname;
    Selector_Log * m_sel_log;

    virtual void Init(const Vec4D *) { }
  public :
    Selector_Base() { m_sel_log=0; m_name=std::string("noname"); }
    virtual ~Selector_Base();
    virtual bool Trigger(const Vec4D *) = 0;
    virtual void Add(Selector_Base *);
    virtual void BuildCuts(Cut_Data *);
    virtual void UpdateCuts(double,double,Cut_Data *);
    virtual void SetRange(std::vector<Flavour>,double,double);
    virtual void SetRange(std::vector<Flavour>,int,double,double);
    virtual void SetRange(std::vector<Flavour>,
			  std::vector<std::pair<double,double> > &);
    virtual void   Output();
    
    virtual bool GetValue(const std::string &name,double &value);

    virtual int    NeedUpdate();
    virtual int    IsConditional();
    virtual void   SetSRange(double _smin,double _smax);
    virtual void   SetName(std::string _name);    
    virtual double ActualValue() const;
    //virtual double Smin();
    //virtual double Smax(); 
    virtual std::string Name(); 
    virtual void SetProcessName(const std::string &name);
    virtual std::string ProcessName() const;
    
    // inline functions
    inline void SetNIn(const size_t nin)   { m_nin=nin;   }
    inline void SetNOut(const size_t nout) { m_nout=nout; }
    inline void SetNTot(const size_t n)    { m_n=n;       }

    inline size_t NIn() const  { return m_nin;  }
    inline size_t NOut() const { return m_nout; }
    inline size_t NTot() const { return m_n;    }

    virtual void AddOnshellCondition(std::string,double) {}
  };

  class Selector_Log {
    std::string   m_name;
    long int      m_rejected, m_passed;
  public :
    Selector_Log(std::string _name) : m_name(_name), m_rejected(0), m_passed(0) { }
    int  Hit();
    int  Hit(bool);
    void Output();
    int  Rejections();
    int  Passed();
  };

  inline int Selector_Log::Hit() { ++m_rejected; return 0; } 
  inline int Selector_Log::Hit(bool hit) { 
    if (hit) { ++m_rejected; return 1; } 
        else { ++m_passed; return 0;}
  }
  inline int  Selector_Log::Rejections() { return m_rejected; }
  inline int  Selector_Log::Passed()     { return m_passed; }
}






/*!
  \class Selector_Base 
  \brief The base class for all selectors.

  This is the base class for all selectors to be built.
  Some of them, the Standard_Selector and the Jet_Finder are
  listed below, as well as the Combined_Selector that is essentially
  a list of selectors.
*/

/*!
  \var   Selector_Base::Selector_Log * sel_log
  \brief The log of the actual selector - keeps track of hits and misses.

  This keeps track of the performance of the selector during its calls, in other 
  words, it is meant to check the performance of integration.
*/  

/*!
  \fn    virtual void Selector_Base::Init(const Vec4D *)
  \brief Boost the vectors in the coorect frame for the Selector to act on them.

  Before cuts are applied in Trigger the system might have to be boosted to 
  a specific frame. This is done here.
*/

/*!
  \fn    virtual bool Selector_Base::Trigger(const Vec4D *)
  \brief The selectors central method - do the vectors pass the selector or not.

  This is the triggering method of the selector. Depending on whether the set of 
  four vectors passes the cuts on the phase space and similar slections or not, the 
  trigger returns 1 or 0 (true or false).
*/

/*!
  \fn    virtual void Selector_Base::Add(Selector_Base *)
  \brief Adding individual selectors to a Combined_Selector

  It is used only by the Combined_Selector to fill the list of individual 
  selectors.
*/

/*!
  \fn    virtual void Selector_Base::BuildCuts(Cut_Data *) 
  \brief Construct a-priori cuts on the phase space.

  This method tries to create a priori cuts on the phase space to be
  used by the phase space integrators. Basically it is tried here
  to minimise the phase space that is maximally allowed for a specific 
  process or group of processes under the cuts applied on it under
  the rule of the slector. It is mandatory to improve the performance
  of the integrator.
*/

/*!
  \fn    virtual void Selector_Base::Output()
  \brief performance checks of the selector.

  This outputs the logs of the individual selectors.
*/




/*!
  \class Selector_Log
  This class is the logbook for the performance of the selector.
  It just counts, how often a given selector returned a zero and
  how often it returned a one.
*/

/*
  \fn    int Selector_Log::Hit(bool hit);
  \brief returns hit and increments rejected or passed counters. 
*/

#endif
