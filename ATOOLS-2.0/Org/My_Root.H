#ifdef ROOT_SUPPORT

#ifndef My_Root_H
#define My_Root_H

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include <map>

#define USING__ATOOLS_Exception

#ifndef USING__My_Root_only
#include "Data_Reader.H"
#include "Exception.H"
#endif

namespace MYROOT {
  
#ifndef USING__My_Root_only
  class My_Root: public ATOOLS::Terminator_Object, public ATOOLS::File_IO_Base {
#else
  class My_Root {
#endif
  public:
      
    typedef std::map<const std::string,TObject *const> String_Object_Map;

  private: 

    String_Object_Map m_objects;

    TApplication *p_root;
    
    TFile *p_file;

    bool m_draw;

  public:

    // constructor
    My_Root(const int argc,char **const argv);

    // destructor
    ~My_Root();

    // member functions
    bool AddObject(TObject *const object,const std::string &key); 
    void Draw();

#ifndef USING__My_Root_only
    void PrepareTerminate();
#endif

    // inline functions
    inline TObject *const GetObject(const std::string &key) const 
    { return m_objects.find(key)->second; }

    inline void operator()(TObject *const object,const std::string &key) 
    { AddObject(object,key); }

    inline TObject *const operator[](const std::string &key) 
    { return GetObject(key); }

    inline void SetDraw(const bool draw) { m_draw=draw; }

    inline bool Draw() const { return m_draw; }

  };// end of class My_Root

  extern My_Root *myroot;
  
}// end of namespace MYROOT

#endif

#endif
