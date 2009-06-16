#ifndef AMEGIC_String_String_Library_H
#define AMEGIC_String_String_Library_H

#include "AMEGIC++/String/MyString.H"

namespace AMEGIC {
  class String_Library {
    static const int buffersize = 512*32;
    int m_mode;
    //std::string pathID,pID;
    void Replace(std::string&,const std::string&,const std::string&);
    int  SingleReplace(std::string&,const std::string&,const std::string&,int&);

    void AddToMakefileAM(std::string, std::string, std::string);
  public:
    String_Library(int mode);
    //String_Library(int,int,Cut_Data*);
    void InitMakefile(std::string);
    void UpdateConfigure(std::string);
    void CreateExtraFiles(std::string);
    int  IsFile(std::string &);
    int  Search(std::string,std::string);
    void Copy(std::string,std::string);
    void AddToMakefile(std::string,std::string,std::string);
  };
}
#endif
