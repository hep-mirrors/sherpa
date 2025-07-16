#ifndef AMEGIC_String_String_Library_H
#define AMEGIC_String_String_Library_H

#include "AMEGIC++/String/MyString.H"

namespace AMEGIC {
  class String_Library {
    static const int buffersize = 512*32;
    int m_mode;
    void Replace(std::string&,const std::string&,const std::string&);
    int  SingleReplace(std::string&,const std::string&,const std::string&,int&);
  public:
    String_Library(int mode);
    void InitMakefile(std::string);
    void UpdateConfigure(std::string);
    int  IsFile(std::string &);
    void Copy(std::string,std::string);
    void AddToCMakefile(std::string,std::string,std::string);
  };
}
#endif
