#ifndef Run_Parameter_H
#define Run_Parameter_H

#include <string>
#include <map>
#ifndef USING__ATOOLS_only
#include "Model_Base.H"
#endif
#include "Data_Read.H"
#include "Flavour.H"
#include "Matrix.H"
#include "Function_Base.H"
#include "MyTiming.H"

#define DEVELOP
#ifndef DEVELOP
#define SAVE_RETURN(x,name) return x;
#else
#include "Data_Return.H"
#define SAVE_RETURN(x,name) return IsDefined(x,name);
#endif

#if defined __GNUC__ && defined __ppc__
#define LIB_SUFFIX ".dylib"
#else
#define LIB_SUFFIX ".so"
#endif

namespace ATOOLS {  

  class Run_Parameter {
    std::string         m_path, m_file;
    static bool s_initialized;
    static std::map<std::string,std::string> s_variables;
  public:
    Run_Parameter();
    ~Run_Parameter();
    void Init(std::string,std::string,int argc=0,char* argv[]=0);
    static void AnalyseEnvironment();
    
    class Gen {  // general parameters
      int                 m_output, m_analysis, m_cutscheme;
      //std::string         m_runtime;
      long int            m_nevents, m_ndicedevents;
      long int            m_seed;
      double              m_ecms, m_ycut, m_delta_r;
      double              m_accu;
      Flavour  m_beam1, m_beam2, p_bunch[2];
      std::string         m_rpa_id, m_flavour_id;
#ifndef USING__ATOOLS_only
      Model_Type::code    m_modelcode;
      MODEL::Model_Base * p_model;
#endif
      MyTiming m_timer;
      double m_timeout;
      bool m_batchmode;
      std::string m_username;
      double m_fac_scale_fac;
      double m_ren_scale_fac;
    public:
      inline long   NumberOfEvents()                      { return m_nevents; }
      inline long   NumberOfDicedEvents()                 { return m_ndicedevents; }
      inline long   RandomSeed()                          { return m_seed;    }
      inline int    Analysis()                            { return m_analysis; }
      inline int    CutScheme()                           { return m_cutscheme; }
      inline void   SetNumberOfEvents(long N)             { m_nevents    = N; }
      inline double Ecms()                                { return m_ecms; }
      void   SetEcms(double _ecms);
      inline double Ycut()                                { return m_ycut; }
      inline void   SetYcut(double _ycut)                 { m_ycut    = _ycut;}
      inline double DeltaR()                              { return m_delta_r; }
      inline void   SetDeltaR(double dr)                  { m_delta_r   = dr;}
      inline const  Flavour Beam1()           { return m_beam1; }
      inline const  Flavour Beam2()           { return m_beam2; }
      void   SetBeam1(const Flavour b);
      void   SetBeam2(const Flavour b);
      inline const  Flavour Bunch(const size_t i)            { return p_bunch[i]; }
      inline void   SetBunch(const Flavour b,const size_t i) { p_bunch[i] = b;    }
      inline void      SetTimeOut(const long int timeout) { m_timeout=timeout;  }
      inline double    TimeOut()                          { return m_timeout;   }
      inline MyTiming &Timer()                            { return m_timer;     }
      inline bool   BatchMode()                           { return m_batchmode; }
      bool CheckTime(const double limit=0.);
      inline void   SetNumberOfDicedEvents(const long diced) { m_ndicedevents=diced; }
      inline const std::string &UserName() { return m_username; }
      inline const std::string &Variable(const std::string &key) { return s_variables[key]; }
#ifndef USING__ATOOLS_only
      void SetScaleFactors(double f, double r) {
	m_fac_scale_fac = f;
	m_ren_scale_fac = r;
      }
      inline double FactorizationScaleFactor() { return m_fac_scale_fac; }
      inline double RenormalizationScaleFactor() { return m_ren_scale_fac; }
      inline void SetModel(MODEL::Model_Base * _mo)       { p_model     = _mo; }
      inline void SetModelType(Model_Type::code _code)    { m_modelcode = _code; }
      inline Model_Type::code Model()                     { return m_modelcode; }

      inline std::string ModelName()                      { return p_model->Name(); } 
      inline int ScalarNumber(const std::string _name)                { return p_model->ScalarNumber(_name); }
      inline double ScalarConstant(const std::string _name)           { return p_model->ScalarConstant(_name); }
      inline double ScalarFunction(const std::string _name,double _t) { return p_model->ScalarFunction(_name,_t); }
      inline double ScalarFunction(const std::string _name)           { return p_model->ScalarFunction(_name); }
      inline CMatrix ComplexMatrix(const std::string _name) { return p_model->ComplexMatrix(_name); }
      inline Complex ComplexMatrixElement(const std::string _name,const int _i,const int _j) { 
	return p_model->ComplexMatrixElement(_name,_i,_j); 
      }
      inline Function_Base * GetScalarFunction(const std::string _name) { 
	return p_model->GetScalarFunction(_name); 
      }
#endif

      /*
      inline int Output()                                 { return m_output; }
      inline int Error()                                  { return (m_output>=0); }
      inline int Events()                                 { return (m_output&1); }
      inline int Tracking()                               { return (m_output&4); }
      inline int Debugging()                              { return (m_output&8); }
      inline int SetOutput(int o)                         { return m_output = o; }
      inline void OutputOff()                             { m_output = -1; }
      */

      inline double Accu()                                { return m_accu; }
      inline std::string GetRpaID()                       { return m_rpa_id; }
      inline std::string GetFlavourID()                   { return m_flavour_id; }
      inline void  SetRpaID(const std::string & _id)      { m_rpa_id = _id; }
      inline void  SetFlavourID(const std::string & _id)  { m_flavour_id = _id; }
      
      friend class Run_Parameter;
    } gen;

    double Picobarn()                       { return 3.89379656e8;} // pb/GeV^2
    double c()                              { return 299792458.e3;} // mm/s
    double hBar()                           { return 6.58211889e-25;} // GeV s
    std::string GetPath()                   { return m_path;} 
    std::string SetPath(std::string path)   { return m_path = path;} 
  };

  extern Run_Parameter rpa;


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
