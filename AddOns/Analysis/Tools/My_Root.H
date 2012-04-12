#ifndef ATOOLS_Org_My_Root_H
#define ATOOLS_Org_My_Root_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__ROOT

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include <map>

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

namespace MYROOT {
  
  class My_Root: public ATOOLS::Terminator_Object, public ATOOLS::File_IO_Base {
  public:
      
    typedef std::map<const std::string,TObject *const> String_Object_Map;

  private: 

    String_Object_Map m_objects;

    TApplication *p_root;
    
    TFile *p_file;

  public:

    // constructor
    My_Root();

    // destructor
    ~My_Root();

    // member functions
    bool AddObject(TObject *const object,const std::string &key); 

    TObject *GetObject(const std::string &key);

    void PrepareTerminate();
    void InitFile();

    std::string GetKey();

    // inline functions
    inline void operator()(TObject *const object,const std::string &key) 
    { AddObject(object,key); }

    inline TObject *operator[](const std::string &key) 
    { return GetObject(key); }

  };// end of class My_Root

  extern My_Root *myroot;
  
}// end of namespace MYROOT

#endif

#endif
