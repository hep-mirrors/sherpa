#ifndef Selector_Base_h
#define Selector_Base_h

#include "Vector.H"
#include "Flavour.H"
#include "Message.H"
#include "Cut_Data.H"

#include <vector>

namespace APHYTOOLS {
  struct Mom_Data {
    std::vector<APHYTOOLS::Flavour> flavs;
    double min,max;
    int    help; 
    int    type;
  };

  class Selector_Data {
  protected :
    std::vector<Mom_Data> data;
  private :
    bool ReadInData(std::string);
  public :
    Selector_Data(std::string path);
    ~Selector_Data()     {};
    int  Size()          { return data.size(); }
    void Data(int,int &,std::vector<APHYTOOLS::Flavour> &,int &,double &,double &);
    void ControlOutput();
  };


  class Selector_Log;

  class Selector_Base {
  protected :
    int            m_n, m_nin, m_nout;
    // long int       m_count;
    Flavour      * m_fl;
    double         m_smin,m_smax;
    std::string    m_name;
    Selector_Log * m_sel_log;

    virtual void Init(const AMATOOLS::Vec4D *) { }
  public :
    Selector_Base() { m_sel_log=0; m_name=std::string("noname"); }
    virtual ~Selector_Base();
    virtual bool Trigger(const AMATOOLS::Vec4D *) = 0;
    virtual void Add(Selector_Base *) {
      AORGTOOLS::msg.Error()<<"Selector_Base::Add : Virtual method."<<std::endl;
    }
    virtual double * ActualValue() {
      AORGTOOLS::msg.Error()<<"Selector_Base::ActualValue :"
			    <<m_name<<" Virtual method."<<std::endl;
    }
    virtual void BuildCuts(Cut_Data *) { 
      AORGTOOLS::msg.Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
    }
    virtual void UpdateCuts(double,double,Cut_Data *) { 
      AORGTOOLS::msg.Error()<<"Selector_Base::BuildCuts : Virtual method."<<std::endl;
    }
    virtual void   SetRange(std::vector<APHYTOOLS::Flavour>,double,double) { 
      AORGTOOLS::msg.Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
    }
    virtual void   SetRange(std::vector<APHYTOOLS::Flavour>,int,double,double) { 
      AORGTOOLS::msg.Error()<<"Selector_Base::SetRange : Virtual method."<<std::endl;
    }
    virtual void   Output();

    virtual void   SetSRange(double _smin,double _smax) { m_smin = _smin; m_smax = _smax; }
    virtual void   SetName(std::string _name)           { m_name = _name; }
    virtual std::string Name()                          { return m_name; }
  };

  class Selector_Log {
    long int passed,rejected;
    std::string   m_name;
  public :
    Selector_Log(std::string _name) : m_name(_name), rejected(0), passed(0) { }
    int  Hit();
    int  Hit(bool);
    void Output();
    int  Rejections();
    int  Passed();
  };

  inline int Selector_Log::Hit() { ++rejected; return 0; } 
  inline int Selector_Log::Hit(bool hit) { 
    if (hit) { ++rejected; return 1; } 
        else { ++passed; return 0;}
  }
  inline int  Selector_Log::Rejections() { return rejected; }
  inline int  Selector_Log::Passed()     { return passed; }
  inline void Selector_Log::Output() { 
    AORGTOOLS::msg.Out()<<"  Selector "<<m_name<<" rejection quota  : "
			<<double(rejected)/double(rejected+passed)
			<<"  ("<<rejected<<" / "<<passed<<")"<<std::endl;
  }
}






/*!
  \class APHYTOOLS::Selector_Base 
  \brief The base class for all selectors.

  This is the base class for all selectors to be built.
  Some of them, the Standard_Selector and the Jet_Finder are
  listed below, as well as the Combined_Selector that is essentially
  a list of selectors.
*/

/*!
  \var   APHYTOOLS::Selector_Base::Selector_Log * sel_log
  \brief The log of the actual selector - keeps track of hits and misses.

  This keeps track of the performance of the selector during its calls, in other 
  words, it is meant to check the performance of integration.
*/  

/*!
  \fn    virtual void APHYTOOLS::Selector_Base::Init(const AMATOOLS::Vec4D *)
  \brief Boost the vectors in the coorect frame for the Selector to act on them.

  Before cuts are applied in Trigger the system might have to be boosted to 
  a specific frame. This is done here.
*/

/*!
  \fn    virtual bool APHYTOOLS::Selector_Base::Trigger(const AMATOOLS::Vec4D *)
  \brief The selectors central method - do the vectors pass the selector or not.

  This is the triggering method of the selector. Depending on whether the set of 
  four vectors passes the cuts on the phase space and similar slections or not, the 
  trigger returns 1 or 0 (true or false).
*/

/*!
  \fn    virtual void APHYTOOLS::Selector_Base::Add(Selector_Base *)
  \brief Adding individual selectors to a Combined_Selector

  It is used only by the Combined_Selector to fill the list of individual 
  selectors.
*/

/*!
  \fn    virtual void APHYTOOLS::Selector_Base::BuildCuts(Cut_Data *) 
  \brief Construct a-priori cuts on the phase space.

  This method tries to create a priori cuts on the phase space to be
  used by the phase space integrators. Basically it is tried here
  to minimise the phase space that is maximally allowed for a specific 
  process or group of processes under the cuts applied on it under
  the rule of the slector. It is mandatory to improve the performance
  of the integrator.
*/

/*!
  \fn    virtual void APHYTOOLS::Selector_Base::Output()
  \brief performance checks of the selector.

  This outputs the logs of the individual selectors.
*/




/*!
  \class APHYTOOLS::Selector_Log
  This class is the logbook for the performance of the selector.
  It just counts, how often a given selector returned a zero and
  how often it returned a one.
*/

/*
  \fn    int APHYTOOLS::Selector_Log::Hit(bool hit);
  \brief returns hit and increments rejected or passed counters. 
*/

#endif
