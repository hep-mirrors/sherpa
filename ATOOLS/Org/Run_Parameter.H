#ifndef ATOOLS_Org_Run_Parameter_H
#define ATOOLS_Org_Run_Parameter_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include <string>
#include <map>
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/MyTiming.H"

namespace ATOOLS {  

  class Run_Parameter {
    std::string         m_path, m_file;
    static bool s_initialized;
    static size_t s_clevel;
    static std::map<std::string,std::string> s_variables;
    static std::vector<std::string> s_cites;
  public:
    Run_Parameter();
    ~Run_Parameter();
    void Init(std::string,std::string,int argc=0,char* argv[]=0);
    static void AnalyseEnvironment();
    
    class Gen {  // general parameters
      int                 m_analysis, m_cutscheme;
      //std::string         m_runtime;
      long int            m_nevents, m_ndicedevents;
      long int            m_seed;
      long int            m_seed2;
      double              m_ecms;
      double              m_cplscale;
      double              m_accu;
      Flavour  m_beam1, m_beam2, p_bunch[2];
      Vec4D m_pbeam[2];
      MyTiming m_timer;
      double m_timeout;
      int m_batchmode;
      std::string m_username;
      bool m_spincorrelations;
    public:
      inline long   NumberOfEvents()                      { return m_nevents; }
      inline long   NumberOfDicedEvents()                 { return m_ndicedevents; }
      inline long   RandomSeed()                          { return m_seed;    }
      inline int    Analysis()                            { return m_analysis; }
      inline int    CutScheme()                           { return m_cutscheme; }
      inline void   SetNumberOfEvents(long N)             { m_nevents    = N; }
      inline double Ecms()                                { return m_ecms; }
      inline double CplScale()                            { return m_cplscale; }
      void   SetEcms(double _ecms);
      void   SetCplScale(double _cplscale);
      inline Vec4D PBeam(short unsigned int i)            { return m_pbeam[i]; }
      void   SetPBeam(short unsigned int i,Vec4D pbeam);
      inline bool   SpinCorrelation()                     { return m_spincorrelations; }
      inline void   SetSpinCorrelation(bool sc)           { m_spincorrelations   = sc;}
      inline const  Flavour Beam1()           { return m_beam1; }
      inline const  Flavour Beam2()           { return m_beam2; }
      void   SetBeam1(const Flavour b);
      void   SetBeam2(const Flavour b);
      inline const  Flavour Bunch(const size_t i)            { return p_bunch[i]; }
      inline void   SetBunch(const Flavour b,const size_t i) { p_bunch[i] = b;    }
      inline void      SetTimeOut(const long int timeout) { m_timeout=timeout;  }
      inline double    TimeOut()                          { return m_timeout;   }
      inline MyTiming &Timer()                            { return m_timer;     }
      inline int    BatchMode()                           { return m_batchmode; }
      bool CheckTime(const double limit=0.);
      inline void   SetNumberOfDicedEvents(const long diced) { m_ndicedevents=diced; }
      inline const std::string &UserName() { return m_username; }
      std::string Variable(const std::string &key,const std::string &def="");
      inline void SetVariable(const std::string &key,const std::string &value) 
      { s_variables[key]=value; }
      inline static const std::vector<std::string> &Citations() { return s_cites; }
      static void AddCitation(const size_t &level,const std::string &cite);

      inline double Accu()                                { return m_accu; }
      
      friend class Run_Parameter;
    } gen;

    double Picobarn()                       { return 3.89379656e8;} // pb/GeV^2
    double c()                              { return 299792458.e3;} // mm/s
    double hBar()                           { return 6.58211889e-25;} // GeV s
    std::string GetPath()                   { return m_path;} 
    std::string SetPath(std::string path)   { return m_path = path;} 
  };

  extern Run_Parameter rpa;

  std::ostream &operator<<(std::ostream &str,const Run_Parameter &rp);


/*!
  \class Run_Parameter
  \brief reads in parameter files and provides the main parameters for the rest of the program

  This class reads in the parameter file Run.dat and provides the program with all
  parameters an swiches necessary for one run.
*/

/*!
  \intern 
  double Picobarn()                       { return 3.89379656e8;} // pb/GeV^2
  double c()                              { return 299792458.e3;} // mm/s
  double hBar()                           { return 6.58211889e-25;} // GeV s
  
  std::string GetPath()                   { return path;} 
  std::string SetPath(std::string _path)  { path = _path;} 
*/

}  // namespace ATOOLS


/*!
 \file
 \brief   contains Run_Parameter and the description the whole ATOOLS package

*/

/*!
  \mainpage

  The ATOOLS package is the basis of the simulation of physics within the SHERPA framework.
  It provides a collection of basic classes and function within the namespace ATOOLS.
  Three main libraries are defined.

  \section amatools Mathematical Tools
  In the directory Math all classes with a pure mathematical meaning are collected.

  \section aphytools Physical Tools
  In the directory Phy all basic classes needed for physics descriptions are collected.
  Here you can find Flavour, Blob, etc.

  \section aorgtools Organisation Tools
  In the directory Org all basic classes needed organisational purposes are collected.
  Here you can find classes for reading in and writing out variables (IO_Handler), accessing
  runtime parameters (Run_Parameter), etc.
*/


/*!
 \namespace ATOOLS
  In the namespace ATOOLS all classes with a pure mathematical meaning are collected.
  In addition all basic classes needed for physics descriptions are collected here, 
  e.g. Flavour, Blob, etc.
  Classes needed for organisational purposes can be found as well, e.g. for reading in and 
  writing out variables (IO_Handler), accessing runtime parameters (Run_Parameter), etc.
*/

#endif  // Run_Parameter_h
