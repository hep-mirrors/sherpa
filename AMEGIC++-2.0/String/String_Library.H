#ifndef String_Library_h
#define String_Library_h

#include "MyString.H"

namespace AMEGIC {
  class String_Library {
    static const int buffersize = 512*32;
    //std::string pathID,pID;
    void Replace(std::string&,const std::string&,const std::string&);
    int  SingleReplace(std::string&,const std::string&,const std::string&,int&);
  public:
    //String_Library(int,int,Cut_Data*);
    void InitMakefile(std::string);
    int  IsFile(std::string &);
    int  Search(std::string,std::string);
    void Copy(std::string,std::string);
    void AddToMakefile(std::string,std::string,std::string);
  };
}
#endif
